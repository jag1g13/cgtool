//
// Created by james on 10/07/15.
//

#ifndef CGTOOL_LIGHT_ARRAY_H
#define CGTOOL_LIGHT_ARRAY_H

template <typename T> class LightArray{
protected:
    T *array_;
    int size_[2] = {0, 0};
    int length_ = 0;

public:
    LightArray<T>(){};
    LightArray<T>(const int x, const int y=1){alloc(x, y);};
    ~LightArray<T>(){
        delete array_;
    };

    LightArray<T>(const LightArray &other);
    LightArray<T>& operator=(const LightArray &other);

    void alloc(const int x, const int y=1);

    T& operator()(const int x, const int y=0);

    const T& at(const int x, const int y=0);

    void smooth(const int n_iter=1);

    void print(const char *format);
};

#endif //CGTOOL_LIGHT_ARRAY_H
