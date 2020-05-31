#include <iostream>
#include <boost/mpi.hpp>
#include <nlohmann/json.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <vector>
#include <utility>
#include <libpng16/png.h>
#include <boost/filesystem.hpp>


namespace bf = boost::filesystem;
using json = nlohmann::json;


template<class T>
class Grid {

    std::vector<T> data;

    size_t data_index(size_t x, size_t y) { return (x * cols) + y; }

public:

    size_t rows, cols;

    Grid(std::vector<T> &matrix, size_t n, size_t m) : rows(n), cols(m), data(std::move(matrix)) {}

    Grid(Grid &grid2) : data(grid2.data), rows(grid2.rows), cols(grid2.cols) {}

    void set(size_t x, size_t y, T value) { data[data_index(x, y)] = value; };

    void swap(Grid &grid2) {
        data.swap(grid2.data);
        rows = grid2.rows;
        cols = grid2.cols;
    }

    auto row_values(size_t row) {
        std::vector<double> shared_row;
        auto begin = data_index(row, 0);
        auto end = data_index(row, cols);
        for (auto i = begin; i < end; ++i) {
            shared_row.emplace_back(data[i]);
        }
        return shared_row;
    }

    void set_row(const std::vector<double> &buffer, size_t row) {
        size_t pos = 0;
        auto begin = data_index(row, 0);
        auto end = data_index(row, cols - 1);
        for (auto i = begin; i < end; ++i) {
            data[i] = buffer[pos];
            ++pos;
        }
    }

    T operator()(size_t x, size_t y) { return data[data_index(x, y)]; }

};


struct conf_data {
    double a;
    size_t dx;
    size_t dy;

    conf_data(size_t xd, size_t yd, double p, double Cp, double k) {
        dx = xd;
        dy = yd;
        a = k / (p * Cp);
    }
};


void read_matrix(std::ifstream &file, std::vector<double> &grid) {
    std::string text;
    auto const start_pos = file.tellg();
    file.ignore(std::numeric_limits<std::streamsize>::max());
    auto const char_count = file.gcount();
    file.seekg(start_pos);
    text = std::string(char_count, char{});
    file.read(&text[0], text.size());
    boost::char_separator<char> sep(" \n\t\r\f");
    boost::tokenizer<boost::char_separator<char>> tokens(text, sep);
    auto p = tokens.begin();
    while (p != tokens.end()) {
        grid.emplace_back(stod(*p));
        ++p;
    }
}


void
update_grid_cells(Grid<double> &old_grid, Grid<double> &new_grid, size_t &start_row, size_t &end_row, conf_data &conf) {
    auto cols_num = old_grid.cols;
    for (auto i = start_row; i < end_row; ++i) {
        for (size_t j = 1; j < cols_num - 1; ++j) {
            double updated_cell =
                    old_grid(i, j) + conf.a * (((old_grid(i - 1, j) - 2 * old_grid(i, j) + old_grid(i + 1, j))
                                                / (conf.dx * conf.dx)) +
                                               ((old_grid(i, j - 1) - 2 * old_grid(i, j) + old_grid(i, j + 1)) /
                                                (conf.dy * conf.dy)));

            new_grid.set(i, j, updated_cell);
        }
    }
}


void save_grid_as_image(std::string image_path, int width, int height, std::vector<std::vector<double>> &matrix) {
    FILE *fp = fopen(image_path.c_str(), "wb");
    if (!fp) {
        fclose(fp);
        throw std::runtime_error("FAILED opening" + image_path);
    }
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) {
        fclose(fp);
        throw std::runtime_error("FAILED PNG configuration");
    }
    png_infop info_ptr;
    if (!(info_ptr = png_create_info_struct(png_ptr)) || setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, nullptr);
        fclose(fp);
    }
    std::vector<uint8_t> image_row(width * height * 3);
    std::vector<uint8_t *> image_row_ptr(height);
    for (int i = 0; i < height; ++i) {
        image_row_ptr[i] = &image_row[i * width * 3];
    }
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    size_t index = 0;
    double max = -1000000;
    for (auto &v: matrix) max = std::max(*std::max_element(v.begin(), v.end()), max);
    double min = 10000000;
    for (auto &v: matrix) min = std::min(*std::min_element(v.begin(), v.end()), min);
    for (const auto &v: matrix) {
        for (const auto &j: v) {
            size_t pos = index * 3;
            double x = 2 * (j - min) / (max - min);
            image_row[pos + 2] = std::max(0.0, 255.0 * (1 - x));
            image_row[pos] = std::max(0.0, 255.0 * (x - 1));
            image_row[pos + 1] = 255 - image_row[pos] - image_row[pos + 2];
            ++index;
        }
    }
    png_set_rows(png_ptr, info_ptr, image_row_ptr.data());
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, nullptr);
    png_write_end(png_ptr, info_ptr);

}


bool von_neumann_criteria(conf_data &conf, size_t dt) {
    return dt > (pow(std::max(conf.dx, conf.dy), 2) / (4 * conf.a));
}


int main(int argc, char *argv[]) {
    std::ifstream file("./../config.json");
    json config_data;
    if (file.is_open()) config_data = json::parse(file);
    size_t rows_num = config_data["n"], cols_num = config_data["m"], dt = config_data["dt"], iter_max = config_data["max_iter_num"];
    std::string save_img_folder = config_data["save_location"];
    conf_data conf(config_data["dx"], config_data["dy"], config_data["p"], config_data["Cp"], config_data["k"]);
    if (von_neumann_criteria(conf, dt)) {
        throw std::runtime_error("von neumann criteria failed");
    }
    boost::mpi::environment env{argc, argv};
    boost::mpi::communicator world;
    auto comm_size = world.size(), rank = world.rank();
    size_t iter_num = 0, when_to_save = config_data["save_iter_img"];
    if (rank != 0) {
        std::ifstream matrix_file(config_data["grid_file_init"]);
        std::vector<double> grid_init;
        read_matrix(matrix_file, grid_init);
        Grid<double> matrix1(grid_init, rows_num, cols_num);
        auto matrix2(matrix1);
        size_t node_row_range = matrix1.rows / (comm_size - 1);
        size_t start_r = (rank - 1) * node_row_range;
        size_t end_r = rank != comm_size - 1 ? rank * node_row_range : matrix1.cols;
        size_t start_rs = rank == 1 ? start_r + 1 : start_r;
        size_t end_rs = rank == comm_size - 1 ? end_r - 1 : end_r;
        while (iter_num < iter_max) {
            update_grid_cells(matrix1, matrix2, start_rs, end_rs, conf);
            if (iter_num % when_to_save == 0) {
                std::vector<double> send_data;
                for (size_t i = start_r; i < end_r; ++i) {
                    for (size_t j = 0; j < matrix1.cols; ++j) {
                        send_data.emplace_back(matrix2(i, j));
                    }
                }
                boost::mpi::gather(world, send_data, 0);
            }
            std::vector<boost::mpi::request> requests;
            if (rank == 1) {
                requests.push_back(world.isend(rank + 1, 123, matrix2.row_values(end_r - 1)));
                std::vector<double> lower_shared;
                requests.push_back(world.irecv(rank + 1, 123, lower_shared));
                boost::mpi::wait_all(requests.begin(), requests.end());
                matrix2.set_row(lower_shared, end_r);
            } else if (rank != comm_size - 1) {
                requests.push_back(world.isend(rank - 1, 123, matrix2.row_values(start_r)));
                requests.push_back(world.isend(rank + 1, 123, matrix2.row_values(end_r - 1)));
                std::vector<double> upper_shared;
                requests.push_back(world.irecv(rank - 1, 123, upper_shared));
                std::vector<double> lower_shared;
                requests.push_back(world.irecv(rank + 1, 123, lower_shared));
                boost::mpi::wait_all(requests.begin(), requests.end());
                matrix2.set_row(upper_shared, start_r - 1);
                matrix2.set_row(lower_shared, end_r);
            } else {
                requests.push_back(world.isend(rank - 1, 123, matrix2.row_values(start_r)));
                std::vector<double> upper_shared;
                requests.push_back(world.irecv(rank - 1, 123, upper_shared));
                boost::mpi::wait_all(requests.begin(), requests.end());
                matrix2.set_row(upper_shared, start_r - 1);
            }
            matrix1.swap(matrix2);
            ++iter_num;
        }
    } else {
        bf::create_directory(bf::path(save_img_folder));
        while (iter_num < iter_max) {
            if (iter_num % when_to_save == 0) {
                std::vector<std::vector<double>> matrix;
                boost::mpi::gather<std::vector<double>>(world, std::vector<double>(), matrix, 0);
                matrix.erase(matrix.begin(), matrix.begin() + 1);
                save_grid_as_image(save_img_folder + "grid" + std::to_string(iter_num) + ".png", cols_num,
                                   rows_num, matrix);
            }
            ++iter_num;
        }
    }
    return 0;
}
