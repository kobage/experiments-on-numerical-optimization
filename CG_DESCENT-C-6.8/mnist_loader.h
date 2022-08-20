#pragma once
#include <vector>
#include <string>
/*
class mnist loader is taken from "https://github.com/arpaka/mnist-loader", under
MIT License

Copyright(c) 2017 Arpaka

Permission is hereby granted, free of charge, to any person obtaining a copy
of this softwareand associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright noticeand this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

class mnist_loader {
private:
    //m_images - - - >  public:
    //std::vector<int> m_labels; - - - > public:
    int m_size;
    int m_rows;
    int m_cols;

    void load_images(std::string file, int num = 0);
    void load_labels(std::string file, int num = 0);
    int  to_int(char* p);

public:
    std::vector<std::vector<double>> m_images;   //from public, to fill X_train and X_test
    std::vector<int> m_labels;   // we madi it public:

    mnist_loader() {}    //------------------ we need it 
    mnist_loader(std::string image_file, std::string label_file, int num);
    mnist_loader(std::string image_file, std::string label_file);
    ~mnist_loader();

    int size() { return m_size; }
    int rows() { return m_rows; }
    int cols() { return m_cols; }

    std::vector<double> images(int id) { return m_images[id]; }
    int labels(int id) { return m_labels[id]; }
};
#include <fstream>
#include <assert.h>

mnist_loader::mnist_loader(std::string image_file,
    std::string label_file,
    int num) :
    m_size(0),
    m_rows(0),
    m_cols(0)
{
    load_images(image_file, num);
    load_labels(label_file, num);
}

mnist_loader::mnist_loader(std::string image_file,
    std::string label_file) :
    mnist_loader(image_file, label_file, 0)
{
    // empty
}

mnist_loader::~mnist_loader()
{
    // empty
}

int
mnist_loader::to_int(char* p)
{
    return ((p[0] & 0xff) << 24) | ((p[1] & 0xff) << 16) |
        ((p[2] & 0xff) << 8) | ((p[3] & 0xff) << 0);
}

void
mnist_loader::load_images(std::string image_file, int num)
{
    std::ifstream ifs(image_file.c_str(), std::ios::in | std::ios::binary);
    char p[4];

    ifs.read(p, 4);
    int magic_number = to_int(p);

    ifs.read(p, 4);
    m_size = to_int(p);
    // limit
    if (num != 0 && num < m_size) m_size = num;

    ifs.read(p, 4);
    m_rows = to_int(p);

    ifs.read(p, 4);
    m_cols = to_int(p);

    char* q = new char[m_rows * m_cols];
    for (int i = 0; i < m_size; ++i) {
        ifs.read(q, m_rows * m_cols);
        std::vector<double> image(m_rows * m_cols);
        for (int j = 0; j < m_rows * m_cols; ++j) {
            image[j] = q[j] / 255.0;
        }
        m_images.push_back(image);
    }
    delete[] q;

    ifs.close();
}

void
mnist_loader::load_labels(std::string label_file, int num)
{
    std::ifstream ifs(label_file.c_str(), std::ios::in | std::ios::binary);
    char p[4];

    ifs.read(p, 4);
    int magic_number = to_int(p);

    ifs.read(p, 4);
    int size = to_int(p);
    // limit
    if (num != 0 && num < m_size) size = num;

    for (int i = 0; i < size; ++i) {
        ifs.read(p, 1);
        int label = p[0];
        m_labels.push_back(label);
    }

    ifs.close();
}