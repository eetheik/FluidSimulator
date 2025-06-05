#include "io.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * @brief A function for printing a density field to the terminal using ANSI color codes.
 *
 * This function takes in a density field @p field with values in [0,1] and renders it by 
 * mapping the field values to one of 216 ANSI color codes, and coloring a space character 
 * with the density associated color. In order for the field to fit in the terminal window it 
 * likely needs to be downsampled. This is handled by the @p step parameter. The field is rendered
 * by looping over the rows and columns, and printing the colored space characters for each entry,
 * and adding a newline character everytime the row ends.
 * 
 * @param n The grid size 
 * @param field A pointer to the @p n * @p n floating point array storing the density field values
 * @param step The step size which is used to downsample the density field so that it fits on the terminal window
 */
void render_colored_field(int n, float *field, int step) {
    // This sets the terminal window start to the top of the screen so that each new frame is rendered
    // so that you don't see the old frame, unless you scroll up ofcourse
    printf("\033[H"); 

    for (int i = 0; i < n; i += step) {
        for (int j = 0; j < n; j += step) {
            // Density values are drawn sideways because the vertical direction has more space at least on my terminal
            // and the size of the space character forces the animation to be a bit stretched in one direction.
            float val = field[i + n * j]; // density values in 0,1 
            val = fmaxf(0.0f, fminf(1.0f, val)); // Clip to 0,1 if not in

            // Map val to an int in 0 - 215, and then add 16 because that's where ANSI colors are stored. 
            // The ANSI colors form a 6x6x6 cube with 216 colors
            int color = 16 + (int)(val * 215);

            // Sets bg color and then draws a space with that color on top of that, 
            // see https://talyian.github.io/ansicolors/
            // \033[0m resets formatting so next element can be printed after this unaffected
            printf("\033[48;5;%dm \033[0m", color); // Print colors
        }
        printf("\n"); // Newline after one finishes
    }

    // Make sure flushes to terminal ASAP
    fflush(stdout);
}


/**
 * @brief A function for initializing a binary file for output.
 *
 * The function initializes a binary file for writing @p nof_frames arrays of size
 * sizeof(float) * @p dim1 * @p dim2 to contiguous disk space. The function takes as
 * input the @p base_filename and adds the time of initialization to the filename.
 * This is so that running a simulation again without changing the filename doesn't override
 * any previous simulation outputs.
 *
 * This function may be useful if one wishes to render the output not to the terminal, but
 * rather using external software, e.g. by generating a GIF in Matplotlib.
 * 
 * @param dim1 An integer denoting the first dimension of the individual arrays being written (we have square arrays so N)
 * @param dim2 An integer denoting the first dimension of the individual arrays being written (we have square arrays so N)
 * @param nof_frames Number of simulation iterations (frames) that are to be saved to the binary output file
 * @param base_filename The path or base filename that the output file will be called. The date is added to the end of this.
 *
 * @return Returns the full initialized filename / filepath, so that it can be passed to other functions e.g. write_array_to_binary
 */
char* init_binary_output_file(const int dim1, const int dim2, const int nof_frames, const char* base_filename) {
    time_t now = time(NULL);
    char timestamp[32];
    char* output_directory = "outputs/";

    // Convert time to string
    strftime(timestamp, sizeof(timestamp), "%d_%m_%Y_%Hh%Mm%Ss", localtime(&now));

    // Assign memory for output_filename 
    size_t output_filename_len = strlen(output_directory) + strlen(base_filename) + strlen(timestamp) + 6;
    char* output_filename = malloc(output_filename_len);
    
    // Create output_filename by concatenating base_filename with timestamp and .bin
    snprintf(output_filename, output_filename_len, "%s%s_%s.bin", output_directory, base_filename, timestamp);

    // Open a binary file in write mode
    FILE* file = fopen(output_filename, "wb");

    // Add an empty byte to the end of the file, forces it to initialize the desired size file
    // This guarantees all writes to the file happen in contiguous disk space. The file isn't
    // growing during the simulation, so no file-handling issues should arise 
    fseek(file, nof_frames * dim1 * dim2 * sizeof(float) - 1, SEEK_SET);
    fputc(0, file);

    fclose(file);
    printf("Simulation output file relative path: %s\n", output_filename);

    return output_filename;
}

/**
 * @brief A function for writing an array to a binary file.
 *
 * The function writes the content of a 2D array (flattened as
 * a 1D array) with dimensions @p dim1 * @p dim2 to a specified binary file.
 * The dimensions are necessary so that the correct sized memory block is
 * written. The passed filename should be the one output by init_binary_output_file.
 *
 * This function may be useful if one wishes to render the output not to the terminal, but
 * rather using external software, e.g. by generating a GIF in Matplotlib.
 *
 * TODO: Prevent misuse, this function could be called on its own without 
 * initialization and that shouldn't work because we open in rb+ mode
 * TODO: Reads and writes are slow so maybe the file could be opened just once, 
 * i.e. initalization opens and this writes? Unsafe?
 *
 * @param array A pointer to the array to be written to the binary file
 * @param current_frame_idx: The current simulation step (frame), this is used to writ to correct part in the file
 * @param dim1 First dimension of the array (number of rows)
 * @param dim2 Second dimension of the array (number of columns)
 * @param filename Path to the output binary (.bin) file obtained from init_binary_output_file
 */
void write_array_to_binary(const float* array, const int current_frame_idx, const int dim1, const int dim2, const char* filename) {
    // "rb+" for reading and writing to binary
    FILE* file = fopen(filename, "rb+");  
    if (!file) {
        fprintf(stderr, "Failed to write to %s\n", filename);
        return;
    }
    fseek(file, current_frame_idx * dim1 * dim2 * sizeof(float), SEEK_SET);
    fwrite(array, sizeof(float), dim1 * dim2, file);
    fclose(file);
}