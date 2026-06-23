package main

import (
	"log"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/gin-gonic/gin"
)

func main() {
	router := gin.Default()
	router.LoadHTMLFiles("template.html")
	router.MaxMultipartMemory = 32 << 20

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
		result, status, errMessage := runEstimate(c)
		if errMessage != "" {
			c.String(status, errMessage)
			return
		}
		c.HTML(http.StatusOK, "template.html", gin.H{
			"Results": result.Results,
		})
	})

	router.POST("/api/estimate", func(c *gin.Context) {
		result, status, errMessage := runEstimate(c)
		if errMessage != "" {
			c.JSON(status, gin.H{"error": errMessage})
			return
		}
		c.JSON(http.StatusOK, result)
	})

	port := os.Getenv("PORT")
	if port == "" {
		port = "8080"
	}
	router.Run(":" + port)
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

func runEstimate(c *gin.Context) (estimateResult, int, string) {
	c.Request.Body = http.MaxBytesReader(c.Writer, c.Request.Body, 64<<20)

	file1, err := c.FormFile("file1")
	if err != nil {
		return estimateResult{}, http.StatusBadRequest, "File1 upload error: " + err.Error()
	}

	file2, err := c.FormFile("file2")
	if err != nil {
		return estimateResult{}, http.StatusBadRequest, "File2 upload error: " + err.Error()
	}

	chain1 := strings.TrimSpace(c.PostForm("chain1"))
	chain2 := strings.TrimSpace(c.PostForm("chain2"))
	bic := c.PostForm("bic")
	exactmo := c.PostForm("exactmo")

	if chain1 == "" || chain2 == "" || bic == "" || exactmo == "" {
		return estimateResult{}, http.StatusBadRequest, "All fields are required"
	}
	if (bic != "bic" && bic != "aic") || (exactmo != "exact" && exactmo != "lh") {
		return estimateResult{}, http.StatusBadRequest, "Invalid estimation options"
	}

	uploadDir, err := os.MkdirTemp("", "bicexact-upload-*")
	if err != nil {
		return estimateResult{}, http.StatusInternalServerError, "Could not prepare upload directory: " + err.Error()
	}
	defer os.RemoveAll(uploadDir)

	path1 := filepath.Join(uploadDir, "protein1"+filepath.Ext(filepath.Base(file1.Filename)))
	path2 := filepath.Join(uploadDir, "protein2"+filepath.Ext(filepath.Base(file2.Filename)))
	if err := c.SaveUploadedFile(file1, path1); err != nil {
		return estimateResult{}, http.StatusInternalServerError, "Could not save file1: " + err.Error()
	}
	if err := c.SaveUploadedFile(file2, path2); err != nil {
		return estimateResult{}, http.StatusInternalServerError, "Could not save file2: " + err.Error()
	}

	log.Printf("File 1 path: %s\n", path1)
	log.Printf("File 2 path: %s\n", path2)
	log.Printf("Chain 1 ID: %s\n", chain1)
	log.Printf("Chain 2 ID: %s\n", chain2)
	log.Printf("BIC: %s\n", bic)
	log.Printf("Exactmo: %s\n", exactmo)

	cmd := exec.Command("./estimate_hinge_numbers", path1, path2, chain1, chain2, bic, exactmo)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return estimateResult{}, http.StatusInternalServerError, "Error executing command: " + err.Error() + "\nOutput: " + string(output)
	}

	return estimateResult{
		Results:   parseResults(string(output)),
		RawOutput: string(output),
	}, http.StatusOK, ""
}

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
