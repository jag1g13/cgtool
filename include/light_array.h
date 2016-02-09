//
// Created by james on 10/07/15.
//

#ifndef CGTOOL_LIGHT_ARRAY_H
#define CGTOOL_LIGHT_ARRAY_H

#include <cassert>
#include <stdexcept>

#include <valarray>
#include <array>

#include "small_functions.h"

template <typename T> class LightArray{
protected:
    std::valarray<T> array_;
    std::array<int, 2> size_ = {{0, 0}};
    const bool safe_;

public:
    LightArray<T>(const bool safe=true) : safe_(safe){};

    LightArray<T>(const int x, const int y = 1, const bool safe = true) : safe_(safe){
        alloc(x, y);
    };

/** \brief Copy constructor */
    LightArray<T>(const LightArray &other) :
            safe_(other.safe_), size_(other.size_), array_(other.array_){
    }

/** \brief Assignment operator */
    LightArray<T> &operator=(const LightArray<T> &other){
        size_ = other.size_;
        // Valarray resizes if necessary since C++11
        array_ = other.array_;
        return *this;
    }

    LightArray<T> &operator/=(const T &div){
        array_ /= div;
        return *this;
    }

    LightArray<T> &operator/=(const LightArray<T> &other){
        assert(size_ == other.size_);
        array_ /= other.array_;
        return *this;
    }

    bool operator==(const LightArray<T> &other) const{
        std::valarray<bool> eq = array_ == other.array_;
        return (size_ == other.size_) && (std::all_of(std::begin(eq), std::end(eq), [](bool b){return b;}));
    }

    void alloc(const int x, const int y = 1){
        size_[0] = x;
        size_[1] = y;
        array_.resize(x * y);
        if(safe_) zero();
    }

    void zero(const T init = 0){
        array_ = init;
    }

    T &operator()(const int x, const int y = 0){
        return array_[x*size_[1] + y];
    }

    const T &at(const int x, const int y = 0) const{
        return array_[(x%size_[0])*size_[1] + (y%size_[1])];
    }

    T sum() const{
        return array_.sum();
    }

    void smooth(const int n_iter=1){
        int jsw = 0;
        int isw = 0;
        for(int ipass = 0; ipass < 2 * n_iter; ipass++){
            jsw = isw;
            for(int i = 1; i < size_[0]-1; i++){
                for(int j = jsw+1; j < size_[1]-1; j += 2){
                    const int loc = i*size_[0] + j;
                    array_[loc] += 0.25 * (array_[loc+size_[0]] + array_[loc+1] + array_[loc-1] + array_[loc-size_[0]] - 4 * array_[loc]);
                }
                jsw = 1 - jsw;
            }
            isw = 1 - isw;
        }
    }

    double mean() const{
        return array_.sum() / array_.size();
    }

/** \brief Print the array - format required to account for different types */
    void print(const char *format="%8.3f") const{
        for(int i = 0; i < size_[0]; i++){
            for(int j = 0; j < size_[1]; j++){
                printf(format, array_[i * size_[0] + j]);
            }
            printf("\n");
        }
    }

    void printCSV(const std::string &filename, const bool suppress_backup=false,
                  const int remove_border=0) const{
        const int r = remove_border;
        assert(r >= 0);
        const std::string file = filename + ".dat";

        // Backup using small_functions.h
        if(!suppress_backup) backup_old_file(file);

        FILE *f = fopen(file.c_str(), "a");
        for(int i=r; i < size_[0]-r; i++){
            for(int j=r; j < size_[1]-r; j++){
                fprintf(f, "%8.3f", array_[i*size_[0] + j]);
            }
            fprintf(f, "\n");
        }

        fclose(f);
    }
};

#endif //CGTOOL_LIGHT_ARRAY_H
