//
// Created by james on 10/07/15.
//

#ifndef CGTOOL_LIGHT_ARRAY_H
#define CGTOOL_LIGHT_ARRAY_H

#include "small_functions.h"

template <typename T> class LightArray{
protected:
    T *array_ = nullptr;
    int size_[2] = {0, 0};
    int length_ = 0;
    bool safe_ = true;

public:
    LightArray<T>(){};
    LightArray<T>(const int x, const int y=1, const bool safe=true){
        alloc(x, y);
        safe_ = safe;
        if(safe_) zero();
    };
    ~LightArray<T>(){
        delete array_;
    };

    LightArray<T>(const LightArray &other);
    LightArray<T>& operator=(const LightArray<T> &other);

    LightArray<T>& operator/=(const int div){
        const double mult = 1./div;
        for(int i=0; i<length_; i++) array_[i] *= mult;
        return *this;
    }

    void alloc(const int x, const int y=1);

    void zero(const T init=0);

    T& operator()(const int x, const int y=0);

    const T& at(const int x, const int y=0) const;

    void smooth(const int n_iter=1);

    double mean() const{
        double mean = 0.;
        for(int i=0; i<length_; i++) mean += array_[i];
        return mean / length_;
    }

    void print(const char *format="%8.3f") const;

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

/** \brief Copy constructor */
template <typename T> LightArray<T>::LightArray(const LightArray &other){
    alloc(other.size_[0], other.size_[1]);
    for(int i=0; i<length_; i++){
        array_[i] = other.array_[i];
    }
}

/** \brief Assignment operator */
template <typename T> LightArray<T>& LightArray<T>::operator=(const LightArray<T> &other){
    if(size_[0]==other.size_[0] && size_[1]==other.size_[1]){
        for(int i=0; i<length_; i++){
            array_[i] = other.array_[i];
        }
    }
}

template <typename T> void LightArray<T>::alloc(const int x, const int y){
    delete array_;
    array_ = nullptr;
    size_[0] = x; size_[1] = y;
    length_ = x * y;
    array_ = new T[length_];
    if(array_ == nullptr) throw std::runtime_error("Could not allocate array");
    zero();
}

template <typename T> void LightArray<T>::zero(const T init){
    for(int i=0; i<length_; i++) array_[i] = init;
}

template <typename T> T& LightArray<T>::operator()(const int x, const int y){
    return array_[x*size_[1] + y];
}

template <typename T> const T& LightArray<T>::at(const int x, const int y) const{
    return array_[(x%size_[0])*size_[1] + (y%size_[1])];
}

/** \brief Apply n iterations of Jacobi smoothing */
template <typename T> void LightArray<T>::smooth(const int n_iter){
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

/** \brief Print the array - format required to account for different types */
template <typename T> void LightArray<T>::print(const char *format) const{
    for(int i=0; i<size_[0]; i++){
        for(int j=0; j<size_[1]; j++){
            printf(format, array_[i*size_[0] + j]);
        }
        printf("\n");
    }
}
#endif //CGTOOL_LIGHT_ARRAY_H
