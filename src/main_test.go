package main

import (
	"bytes"
	"context"
	"encoding/json"
	"io"
	"mime/multipart"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func testService(t *testing.T, script string) *estimateService {
	t.Helper()
	executable := filepath.Join(t.TempDir(), "estimator")
	if err := os.WriteFile(executable, []byte("#!/bin/sh\n"+script), 0700); err != nil {
		t.Fatal(err)
	}
	return &estimateService{slots: make(chan struct{}, 1), timeout: time.Second, executable: executable}
}

func uploadRequest(t *testing.T, route, filename, chain string, fileSize int) *http.Request {
	t.Helper()
	var body bytes.Buffer
	w := multipart.NewWriter(&body)
	for _, name := range []string{"file1", "file2"} {
		part, err := w.CreateFormFile(name, filename)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := io.CopyN(part, strings.NewReader(strings.Repeat("X", fileSize)), int64(fileSize)); err != nil {
			t.Fatal(err)
		}
	}
	for name, value := range map[string]string{"chain1": chain, "chain2": "A", "bic": "bic", "exactmo": "exact"} {
		if err := w.WriteField(name, value); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	req := httptest.NewRequest(http.MethodPost, route, &body)
	req.Header.Set("Content-Type", w.FormDataContentType())
	return req
}

func TestEstimateSuccessAndCleanup(t *testing.T) {
	tmp := t.TempDir()
	t.Setenv("TMPDIR", tmp)
	s := testService(t, "test -f \"$1\" && test -f \"$2\" || exit 1\nprintf 'Protein length: 99\nThe number of hinges is: 1\n'")
	router := newRouter(s)
	for _, route := range []string{"/api/estimate", "/upload"} {
		t.Run(route, func(t *testing.T) {
			req := uploadRequest(t, route, "../../protein.PDB", "A", 2<<20)
			response := httptest.NewRecorder()
			router.ServeHTTP(response, req)
			if response.Code != http.StatusOK || !strings.Contains(response.Body.String(), "99") {
				t.Fatalf("unexpected response: %d %s", response.Code, response.Body.String())
			}
			if route == "/api/estimate" {
				var result estimateResult
				if err := json.Unmarshal(response.Body.Bytes(), &result); err != nil || result.Results["Protein length"] != "99" {
					t.Fatalf("invalid result: %+v, %v", result, err)
				}
			}
			entries, err := os.ReadDir(tmp)
			if err != nil || len(entries) != 0 {
				t.Fatalf("temporary uploads were not removed: %v, %v", entries, err)
			}
		})
	}
}

func TestRejectInvalidUploads(t *testing.T) {
	router := newRouter(testService(t, "exit 99"))
	for _, tc := range []struct {
		name, filename, chain string
		size                  int
	}{
		{"extension", "protein.txt", "A", 1},
		{"empty", "protein.pdb", "A", 0},
		{"empty chain", "protein.pdb", "", 1},
		{"long chain", "protein.pdb", "ABCDEFGHI", 1},
		{"control character", "protein.pdb", "A\nB", 1},
		{"shell characters", "protein.pdb", "A;id", 1},
	} {
		t.Run(tc.name, func(t *testing.T) {
			response := httptest.NewRecorder()
			router.ServeHTTP(response, uploadRequest(t, "/api/estimate", tc.filename, tc.chain, tc.size))
			if response.Code != http.StatusBadRequest {
				t.Fatalf("status = %d: %s", response.Code, response.Body.String())
			}
		})
	}
}

type repeatingReader struct{}

func (repeatingReader) Read(p []byte) (int, error) {
	for i := range p {
		p[i] = 'X'
	}
	return len(p), nil
}

func TestOversizedUpload(t *testing.T) {
	tmp := t.TempDir()
	t.Setenv("TMPDIR", tmp)
	router := newRouter(testService(t, "exit 99"))
	prefix := "--boundary\r\nContent-Disposition: form-data; name=\"file1\"; filename=\"x.pdb\"\r\n\r\n"
	reader := io.MultiReader(strings.NewReader(prefix), io.LimitReader(repeatingReader{}, maxUploadBytes+1))
	req := httptest.NewRequest(http.MethodPost, "/api/estimate", reader)
	req.Header.Set("Content-Type", "multipart/form-data; boundary=boundary")
	response := httptest.NewRecorder()
	router.ServeHTTP(response, req)
	if response.Code != http.StatusRequestEntityTooLarge {
		t.Fatalf("status = %d: %s", response.Code, response.Body.String())
	}
	entries, err := os.ReadDir(tmp)
	if err != nil || len(entries) != 0 {
		t.Fatalf("partial upload was not removed: %v, %v", entries, err)
	}
}

func TestEstimatorFailures(t *testing.T) {
	for _, tc := range []struct {
		name, script string
		status       int
	}{
		{"internal error", "echo /private/server/secret >&2; exit 1", http.StatusInternalServerError},
		{"invalid structure", "echo /private/server/secret >&2; exit 2", http.StatusBadRequest},
		{"timeout", "exec sleep 10", http.StatusGatewayTimeout},
		{"excessive output", "exec yes /private/server/secret", http.StatusInternalServerError},
	} {
		t.Run(tc.name, func(t *testing.T) {
			s := testService(t, tc.script)
			s.timeout = 2 * time.Second
			if tc.name == "timeout" {
				s.timeout = 100 * time.Millisecond
			}
			response := httptest.NewRecorder()
			start := time.Now()
			newRouter(s).ServeHTTP(response, uploadRequest(t, "/api/estimate", "p.pdb", "A", 1))
			if response.Code != tc.status || strings.Contains(response.Body.String(), "/private/") {
				t.Fatalf("unexpected response: %d %s", response.Code, response.Body.String())
			}
			if time.Since(start) > 3*time.Second || len(s.slots) != 0 {
				t.Fatal("estimator did not release resources promptly")
			}
		})
	}
}

func TestBusyServerRejectsBeforeReadingUpload(t *testing.T) {
	s := testService(t, "exit 99")
	s.slots <- struct{}{}
	response := httptest.NewRecorder()
	newRouter(s).ServeHTTP(response, httptest.NewRequest(http.MethodPost, "/api/estimate", nil))
	if response.Code != http.StatusServiceUnavailable || response.Header().Get("Retry-After") == "" {
		t.Fatalf("unexpected response: %d", response.Code)
	}
}

func TestCanceledRequestStopsEstimator(t *testing.T) {
	s := testService(t, "exec sleep 10")
	req := uploadRequest(t, "/api/estimate", "p.pdb", "A", 1)
	ctx, cancel := context.WithCancel(req.Context())
	cancel()
	start := time.Now()
	response := httptest.NewRecorder()
	newRouter(s).ServeHTTP(response, req.WithContext(ctx))
	if response.Code != http.StatusInternalServerError || time.Since(start) > time.Second || len(s.slots) != 0 {
		t.Fatal("canceled request did not release resources")
	}
}

func TestCrossOriginUploadsBlocked(t *testing.T) {
	router := newRouter(testService(t, "exit 99"))
	for _, route := range []string{"/upload", "/api/estimate"} {
		req := uploadRequest(t, route, "p.pdb", "A", 1)
		req.Header.Set("Origin", "https://attacker.example")
		response := httptest.NewRecorder()
		router.ServeHTTP(response, req)
		if response.Code != http.StatusForbidden {
			t.Fatalf("cross-origin status = %d", response.Code)
		}
	}
}
