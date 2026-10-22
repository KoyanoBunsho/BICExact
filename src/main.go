package main

import (
	"bytes"
	"context"
	"errors"
	"io"
	"log"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"time"

	"github.com/gin-gonic/gin"
)

func main() {
	port := os.Getenv("PORT")
	if port == "" {
		port = "8080"
	}
	server := &http.Server{
		Addr: ":" + port,
		Handler: newRouter(&estimateService{
			slots: make(chan struct{}, 2), timeout: 30 * time.Second,
			executable: "./estimate_hinge_numbers",
		}),
		ReadHeaderTimeout: 5 * time.Second,
		ReadTimeout:       15 * time.Second,
		WriteTimeout:      50 * time.Second,
		IdleTimeout:       60 * time.Second,
		MaxHeaderBytes:    16 << 10,
	}
	log.Fatal(server.ListenAndServe())
}

func newRouter(estimator *estimateService) http.Handler {
	router := gin.Default()
	if err := router.SetTrustedProxies(nil); err != nil {
		panic(err)
	}
	router.Use(func(c *gin.Context) {
		c.Header("X-Content-Type-Options", "nosniff")
		c.Header("X-Frame-Options", "DENY")
		c.Header("Referrer-Policy", "no-referrer")
		c.Header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")
		c.Next()
	})
	router.LoadHTMLFiles("template.html")

	distDir := frontendDistDir()
	if _, err := os.Stat(filepath.Join(distDir, "index.html")); err == nil {
		router.Static("/assets", filepath.Join(distDir, "assets"))
		router.NoRoute(func(c *gin.Context) {
			c.File(filepath.Join(distDir, "index.html"))
		})
	} else {
		router.GET("/", func(c *gin.Context) {
			c.HTML(http.StatusOK, "template.html", nil)
		})
	}

	router.POST("/upload", func(c *gin.Context) {
		result, status, errMessage := estimator.runEstimate(c)
		if errMessage != "" {
			c.String(status, errMessage)
			return
		}
		c.HTML(http.StatusOK, "template.html", gin.H{
			"Results": result.Results,
		})
	})

	router.POST("/api/estimate", func(c *gin.Context) {
		result, status, errMessage := estimator.runEstimate(c)
		if errMessage != "" {
			c.JSON(status, gin.H{"error": errMessage})
			return
		}
		c.JSON(http.StatusOK, result)
	})

	return http.NewCrossOriginProtection().Handler(router)
}

type estimateResult struct {
	Results   map[string]string `json:"results"`
	RawOutput string            `json:"rawOutput"`
}

func frontendDistDir() string {
	if dir := os.Getenv("FRONTEND_DIST"); dir != "" {
		return dir
	}
	return filepath.Join("..", "frontend", "dist")
}

const maxUploadBytes = 64 << 20
const maxOutputBytes = 64 << 10

var chainPattern = regexp.MustCompile(`^[A-Za-z0-9]{1,8}$`)

type estimateService struct {
	slots      chan struct{}
	timeout    time.Duration
	executable string
}

func (s *estimateService) runEstimate(c *gin.Context) (estimateResult, int, string) {
	// Reserve capacity before parsing uploads, which may consume memory and disk.
	select {
	case s.slots <- struct{}{}:
		defer func() { <-s.slots }()
	default:
		c.Header("Retry-After", "5")
		return estimateResult{}, http.StatusServiceUnavailable, "Server is busy; please retry shortly"
	}
	c.Request.Body = http.MaxBytesReader(c.Writer, c.Request.Body, maxUploadBytes)
	if err := c.Request.ParseMultipartForm(1 << 20); err != nil {
		var tooLarge *http.MaxBytesError
		if errors.As(err, &tooLarge) {
			return estimateResult{}, http.StatusRequestEntityTooLarge, "Upload exceeds the 64 MiB limit"
		}
		return estimateResult{}, http.StatusBadRequest, "Invalid multipart upload"
	}
	defer c.Request.MultipartForm.RemoveAll()

	file1, err := c.FormFile("file1")
	if err != nil {
		return estimateResult{}, http.StatusBadRequest, "File1 is required"
	}

	file2, err := c.FormFile("file2")
	if err != nil {
		return estimateResult{}, http.StatusBadRequest, "File2 is required"
	}
	ext1, ext2 := strings.ToLower(filepath.Ext(file1.Filename)), strings.ToLower(filepath.Ext(file2.Filename))
	validExtension := func(ext string) bool { return ext == ".pdb" || ext == ".ent" || ext == ".cif" }
	if !validExtension(ext1) || !validExtension(ext2) || file1.Size == 0 || file2.Size == 0 {
		return estimateResult{}, http.StatusBadRequest, "Upload two nonempty .pdb, .ent, or .cif files"
	}

	chain1 := strings.TrimSpace(c.PostForm("chain1"))
	chain2 := strings.TrimSpace(c.PostForm("chain2"))
	bic := c.PostForm("bic")
	exactmo := c.PostForm("exactmo")

	if !chainPattern.MatchString(chain1) || !chainPattern.MatchString(chain2) {
		return estimateResult{}, http.StatusBadRequest, "Chain IDs must contain 1 to 8 ASCII letters or digits"
	}
	if (bic != "bic" && bic != "aic") || (exactmo != "exact" && exactmo != "lh") {
		return estimateResult{}, http.StatusBadRequest, "Invalid estimation options"
	}

	uploadDir, err := os.MkdirTemp("", "bicexact-upload-*")
	if err != nil {
		log.Printf("Prepare upload directory: %v", err)
		return estimateResult{}, http.StatusInternalServerError, "Could not prepare upload directory"
	}
	defer os.RemoveAll(uploadDir)

	path1 := filepath.Join(uploadDir, "protein1"+ext1)
	path2 := filepath.Join(uploadDir, "protein2"+ext2)
	if err := c.SaveUploadedFile(file1, path1); err != nil {
		log.Printf("Save file1: %v", err)
		return estimateResult{}, http.StatusInternalServerError, "Could not save file1"
	}
	if err := c.SaveUploadedFile(file2, path2); err != nil {
		log.Printf("Save file2: %v", err)
		return estimateResult{}, http.StatusInternalServerError, "Could not save file2"
	}

	ctx, cancel := context.WithTimeout(c.Request.Context(), s.timeout)
	defer cancel()
	cmd := exec.CommandContext(ctx, s.executable, path1, path2, chain1, chain2, bic, exactmo)
	cmd.WaitDelay = time.Second
	output := &limitedOutput{cancel: cancel}
	cmd.Stdout, cmd.Stderr = output, output
	if err := cmd.Run(); err != nil {
		if errors.Is(ctx.Err(), context.DeadlineExceeded) {
			return estimateResult{}, http.StatusGatewayTimeout, "Estimation exceeded the time limit"
		}
		var exitErr *exec.ExitError
		if errors.As(err, &exitErr) && exitErr.ExitCode() == 2 {
			return estimateResult{}, http.StatusBadRequest, "Invalid structures: check coordinates, chain IDs, and the 10000-residue limit"
		}
		log.Printf("Estimator failed: %v", err)
		return estimateResult{}, http.StatusInternalServerError, "Estimation failed"
	}

	return estimateResult{
		Results:   parseResults(output.String()),
		RawOutput: output.String(),
	}, http.StatusOK, ""
}

// exec.Cmd serializes writes when Stdout and Stderr share the same writer.
type limitedOutput struct {
	buffer bytes.Buffer
	cancel context.CancelFunc
}

func (b *limitedOutput) Write(p []byte) (int, error) {
	if len(p) > maxOutputBytes-b.buffer.Len() {
		b.cancel()
		return 0, io.ErrShortBuffer
	}
	return b.buffer.Write(p)
}

func (b *limitedOutput) String() string { return b.buffer.String() }

func parseResults(output string) map[string]string {
	results := make(map[string]string)
	for _, line := range strings.Split(output, "\n") {
		parts := strings.SplitN(line, ":", 2)
		if len(parts) == 2 {
			results[strings.TrimSpace(parts[0])] = strings.TrimSpace(parts[1])
		}
	}
	return results
}
