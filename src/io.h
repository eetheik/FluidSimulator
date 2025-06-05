#ifndef IO_H
#define IO_H

void render_colored_field(int n, float *field, int step);
char* init_binary_output_file(const int dim1, const int dim2, const int nof_frames, const char* base_filename);
void write_array_to_binary(const float* array, const int current_frame_idx, const int dim1, const int dim2, const char* filename); 

#endif
