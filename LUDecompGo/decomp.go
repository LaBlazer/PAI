package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"
)

var wg sync.WaitGroup

func loadMatrix(filename string) [][]float64 {
	file, err := os.OpenFile(filename, os.O_RDONLY, os.ModePerm)

	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// get mtx size
	scanner.Scan()
	tokens := strings.Split(scanner.Text(), " ")
	n := len(tokens)

	// allocate mtx
	mtx := make([][]float64, n)
	rows := make([]float64, n*n)
	for i := 0; i < n; i++ {
		mtx[i] = rows[i*n : (i+1)*n]
	}

	// convert text to mtx data
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			mtx[i][j], _ = strconv.ParseFloat(tokens[j], 64)
		}

		if !scanner.Scan() {
			break
		}

		tokens = strings.Split(scanner.Text(), " ")
	}

	if err := scanner.Err(); err != nil {
		log.Fatalf("scan file error: %v", err)
	}

	return mtx
}

func printMatrix(mtx [][]float64) {
	if len(mtx) > 10 {
		return
	}

	for i := 0; i < len(mtx); i++ {
		for j := 0; j < len(mtx); j++ {
			fmt.Print("\t", mtx[i][j])
		}
		fmt.Println()
	}
	fmt.Println()
}

func gaussianSerial(mtx [][]float64) {
	n := len(mtx)

	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			lik := mtx[i][k] / mtx[k][k]

			for j := k + 1; j < n; j++ {
				mtx[i][j] -= lik * mtx[k][j]
			}

			mtx[i][k] = lik
		}
	}

	printMatrix(mtx)

}
func gaussianParallel(mtx [][]float64, threadCount int) {
	n := len(mtx)

	for k := 0; k < n; k++ {
		wg.Add(threadCount)
		for t := 0; t < threadCount; t++ {
			go func(threadId int) {
				for i := k + int(threadId) + 1; i < n; i += threadCount {
					lik := mtx[i][k] / mtx[k][k]

					for j := k + 1; j < n; j++ {
						mtx[i][j] -= lik * mtx[k][j]
					}

					mtx[i][k] = lik
				}

				wg.Done()
			}(t)
		}

		wg.Wait()
	}

	printMatrix(mtx)
}

func main() {
	filename := "micro.txt"
	threads := 4

	fmt.Println("Loading from ", filename)

	mtx := loadMatrix(filename)

	fmt.Println("Input mtx size ", len(mtx))
	printMatrix(mtx)

	fmt.Println("Gaussian serial")

	//gaussianSerial(mtx)

	fmt.Println("Gaussian parallel")

	gaussianParallel(mtx, threads)
}
