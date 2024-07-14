package main

import (
	"log"
	"net/http"
	"os/exec"
	"strings"

	"github.com/gin-gonic/gin"
)

func main() {
	router := gin.Default()
	router.LoadHTMLFiles("template.html")

	router.GET("/", func(c *gin.Context) {
		c.HTML(http.StatusOK, "template.html", nil)
	})

	router.POST("/upload", func(c *gin.Context) {
		file1, err := c.FormFile("file1")
		if err != nil {
			c.String(http.StatusBadRequest, "File1 upload error: %s", err.Error())
			return
		}

		file2, err := c.FormFile("file2")
		if err != nil {
			c.String(http.StatusBadRequest, "File2 upload error: %s", err.Error())
			return
		}

		chain1 := c.PostForm("chain1")
		chain2 := c.PostForm("chain2")
		bic := c.PostForm("bic")
		exactmo := c.PostForm("exactmo")

		if chain1 == "" || chain2 == "" || bic == "" || exactmo == "" {
			c.String(http.StatusBadRequest, "All fields are required")
			return
		}

		path1 := file1.Filename
		path2 := file2.Filename
		if err := c.SaveUploadedFile(file1, path1); err != nil {
			c.String(http.StatusInternalServerError, "Could not save file1: %s", err.Error())
			return
		}
		if err := c.SaveUploadedFile(file2, path2); err != nil {
			c.String(http.StatusInternalServerError, "Could not save file2: %s", err.Error())
			return
		}
		log.Printf("File 1 path: %s\n", path1)
		log.Printf("File 2 path: %s\n", path2)
		log.Printf("Chain 1 ID: %s\n", chain1)
		log.Printf("Chain 2 ID: %s\n", chain2)
		log.Printf("BIC: %s\n", bic)
		log.Printf("Exactmo: %s\n", exactmo)

		cmd := exec.Command("docker", "exec", "mplh", "./estimate_hinge_numbers", path1, path2, chain1, chain2, bic, exactmo)
		output, err := cmd.CombinedOutput()
		if err != nil {
			c.String(http.StatusInternalServerError, "Error executing command: %s\nOutput: %s", err.Error(), string(output))
			return
		}

		outputLines := strings.Split(string(output), "\n")
		results := make(map[string]string)
		for _, line := range outputLines {
			parts := strings.SplitN(line, ":", 2)
			if len(parts) == 2 {
				results[strings.TrimSpace(parts[0])] = strings.TrimSpace(parts[1])
			}
		}

		c.HTML(http.StatusOK, "template.html", gin.H{
			"Results": results,
		})
	})

	router.Run(":8080")
}
