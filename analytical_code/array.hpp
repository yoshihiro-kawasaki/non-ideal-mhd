#ifndef ARRAY_HPP_
#define ARRAY_HPP_

/*
    2023/11/17 kawasaki 
*/

#include <iostream>
#include <cmath>
#include <string>
#include <complex>

namespace array
{
/*
#############################################################
define
#############################################################
*/

// #define CHECK_ARGUMENT_
// #define CHECK_BOUNDS_


/*
#############################################################
enum ArrayStatus
#############################################################
*/

enum class ArrayStatus {
    empty,     // not allocate array
    allocated  // allocate array
};


/*
#############################################################
class Array1D
#############################################################
*/

template<typename T>
class Array1D
{
public:
    // constructor
    Array1D();
    explicit Array1D(const int n1);
    Array1D(const int n1, const T &a);
    Array1D(const int n1, const T *a);
    Array1D(const Array1D &rhs);

    // destructor
    ~Array1D();

    // operator
    Array1D & operator=(const Array1D &rhs);
    Array1D & operator=(const T &a);
    inline T & operator[](const int i);
    inline const T & operator[](const int i) const;

    // menmber function
	inline int GetDim1() const { return n1_;}
    inline int Size() const { return n1_; }

	void Resize(const int n1); 
	void Assign(const int n1, const T &a);

private:
    int n1_;
    T* pv_;
    ArrayStatus status_;

    void AllocateArray();
    void DeleteArray();
};

// constructor

template<typename T>
Array1D<T>::Array1D() : n1_(0), pv_(nullptr), status_(ArrayStatus::empty) { }

template<typename T>
Array1D<T>::Array1D(const int n1)
    : n1_(n1), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
}

template<typename T>
Array1D<T>::Array1D(const int n1, const T &a)
    : n1_(n1), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) pv_[i] = a;
}

template<typename T>
Array1D<T>::Array1D(const int n1, const T *a)
    : n1_(n1), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) pv_[i] = *a++;
}

template<typename T>
Array1D<T>::Array1D(const Array1D<T> &rhs)
    : n1_(rhs.n1_), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1_; ++i) pv_[i] = rhs[i];
}

// destructor

template<typename T>
Array1D<T>::~Array1D()
{
    DeleteArray();
}

// operator

template<typename T>
Array1D<T> & Array1D<T>::operator=(const Array1D<T> &rhs)
{
    if (this != &rhs) {

        if (n1_ != rhs.n1_) {
            DeleteArray();
            n1_ = rhs.n1_;
            AllocateArray();
        }

        for (int i = 0; i < n1_; ++i) pv_[i] = rhs[i];

    }

    return *this;
}


template<typename T>
Array1D<T> & Array1D<T>::operator=(const T &a)
{
    if (status_ == ArrayStatus::allocated) {
        for (int i = 0; i < n1_; ++i) pv_[i] = a;
    }
    return *this;
}

template<typename T>
inline T & Array1D<T>::operator[](const int i)
{
    return pv_[i];
}


template<typename T>
inline const T & Array1D<T>::operator[](const int i) const
{
    return pv_[i];
}

// member function


template<typename T>
void Array1D<T>::Resize(const int n1)
{
    if (n1 != n1_) {
        DeleteArray();
        n1_ = n1;
        AllocateArray();
    }
}

template<typename T>
void Array1D<T>::Assign(const int n1, const T &a)
{
    if (n1 != n1_) {
        DeleteArray();
        n1_ = n1;
        AllocateArray();
    }

    for (int i = 0; i < n1; ++i) pv_[i] = a;
}

template<typename T>
void Array1D<T>::AllocateArray()
{
    if (n1_ > 0 && status_ == ArrayStatus::empty) {
        pv_ = new T[n1_];
        status_ = ArrayStatus::allocated;
    }
}

template<typename T>
void Array1D<T>::DeleteArray()
{
    if (status_ == ArrayStatus::allocated) {
        delete [] pv_;
        pv_ = nullptr;
        status_ = ArrayStatus::empty;
    }
}


/*
#############################################################
class Array2D
#############################################################
*/


template<typename T>
class Array2D
{
public:
    // constructor
    Array2D();
    Array2D(const int n1, const int n2);
    Array2D(const int n1, const int n2, const T &a);
    Array2D(const int n1, const int n2, const T *a);
    Array2D(const Array2D &rhs);

    // destructor
    ~Array2D();

    // operator
    Array2D & operator=(const Array2D &rhs);
    Array2D & operator=(const T &a);

	// typedef T value_type; 
	inline T* operator[](const int i);
	inline const T* operator[](const int i) const;

    // menmber function
	inline int GetDim1() const { return n1_;};
	inline int GetDim2() const { return n2_; };

	void Resize(const int n1, const int n2); 
	void Assign(const int n1, const int n2, const T &a); 

private:
    int n1_;
    int n2_;
    T **pv_;
    ArrayStatus status_;

    void AllocateArray();
    void DeleteArray();
};

// constructor

template<typename T>
Array2D<T>::Array2D() : n1_(0), n2_(0), pv_(nullptr), status_(ArrayStatus::empty) { }

template<typename T>
Array2D<T>::Array2D(const int n1, const int n2)
    : n1_(n1), n2_(n2), pv_(nullptr), status_(ArrayStatus::empty)
{
#ifdef CHECK_ARGUMENT_
    if (n1 < 0) {
        std::cout << "n1 < 0" << std::endl;
    }
    if (n2 < 0) {
        std::cout << "n1 < 0" << std::endl;
    }
#endif
    AllocateArray();
}

template<typename T>
Array2D<T>::Array2D(const int n1, const int n2, const T &a)
    : n1_(n1), n2_(n2), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            pv_[i][j] = a;
        }
    }
}

template<typename T>
Array2D<T>::Array2D(const int n1, const int n2, const T *a)
    : n1_(n1), n2_(n2), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            pv_[i][j] = *a++;
        }
    }
}

template<typename T>
Array2D<T>::Array2D(const Array2D &rhs)
    : n1_(rhs.n1_), n2_(rhs.n2_), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            pv_[i][j] = rhs[i][j];
        }
    }
}

// destructor
template<typename T>
Array2D<T>::~Array2D()
{
    DeleteArray();
}

// operator
template<typename T>
Array2D<T> & Array2D<T>::operator=(const Array2D &rhs)
{
    if (this != &rhs) {

        if (n1_ != rhs.n1_ || n2_ != rhs.n2_) {
            DeleteArray();
            n1_ = rhs.n1_;
            n2_ = rhs.n2_;
            AllocateArray();
        }

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                pv_[i][j] = rhs[i][j];
            }
        }

    }

    return *this;
}

template<typename T>
Array2D<T> & Array2D<T>::operator=(const T &a)
{
    if (status_ == ArrayStatus::allocated) {
        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                pv_[i][j] = a;
            }
        }
    }

    return *this;
}

template<typename T>
inline T* Array2D<T>::operator[](const int i)
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

template<typename T>
inline const T* Array2D<T>::operator[](const int i) const
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// menmber function

template<typename T>
void Array2D<T>::Resize(const int n1, const int n2)
{
    if (n1 != n1_ || n2 != n2_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        AllocateArray();
    }
}

template<typename T>
void Array2D<T>::Assign(const int n1, const int n2, const T &a)
{
    if (n1 != n1_ || n2 != n2_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        AllocateArray();
    }
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            pv_[i][j] = a;
        }
    }
}



template<typename T>
void Array2D<T>::AllocateArray()
{
    if (status_ == ArrayStatus::empty) {
        int size = n1_ * n2_;
        pv_    = new T*[n1_];
        pv_[0] = new T[size];
        for (int i = 0; i < n1_; ++i) {
            pv_[i] = pv_[0] + i*n2_;
        }
        status_ = ArrayStatus::allocated;
    }
}


template<typename T>
void Array2D<T>::DeleteArray()
{
    if (status_ == ArrayStatus::allocated) {
        delete [] pv_[0];
        delete [] pv_;
        pv_ = nullptr;
        status_ = ArrayStatus::empty;
    }
}

/*
#############################################################
class Array3D
#############################################################
*/

template<typename T>
class Array3D
{
public:
    // constructor
    Array3D();
    Array3D(const int n1, const int n2, const int n3);
    Array3D(const int n1, const int n2, const int n3, const T &a);
    Array3D(const int n1, const int n2, const int n4, const T *a);
    Array3D(const Array3D &rhs);

    // destructor
    ~Array3D();

    // operator
    Array3D & operator=(const Array3D &rhs);
    Array3D & operator=(const T &a);
	// typedef T value_type; 
	inline T** operator[](const int i);
	inline const T* const * operator[](const int i) const;

    // menmber function
	inline int GetDim1() const { return n1_; }
	inline int GetDim2() const { return n2_; }
    inline int GetDim3() const { return n3_; }

	void Resize(const int n1, const int n2, const int n3); 
	void Assign(const int n1, const int n2, const int n3, const T &a); 

private:
    int n1_;
    int n2_;
    int n3_;
    T ***pv_;
    ArrayStatus status_;

    void AllocateArray();
    void DeleteArray();
};

// constructor

template<typename T>
Array3D<T>::Array3D() : n1_(0), n2_(0), n3_(0), pv_(nullptr), status_(ArrayStatus::empty) { }

template<typename T>
Array3D<T>::Array3D(const int n1, const int n2, const int n3)
    : n1_(n1), n2_(n2), n3_(n3), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
}

template<typename T>
Array3D<T>::Array3D(const int n1, const int n2, const int n3, const T &a)
    : n1_(n1), n2_(n2), n3_(n3), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                pv_[i][j][k] = a;
            }
        }
    }
}


template<typename T>
Array3D<T>::Array3D(const int n1, const int n2, const int n3, const T *a)
    : n1_(n1), n2_(n2), n3_(n3), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                pv_[i][j][k] = *a++;
            }
        }
    }
}


template<typename T>
Array3D<T>::Array3D(const Array3D<T> &rhs)
    : n1_(rhs.n1_), n2_(rhs.n2_), n3_(rhs.n3_), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                pv_[i][j][k] = rhs[i][j][k];
            }
        }
    }
}

// destructor
template<typename T>
Array3D<T>::~Array3D()
{
    DeleteArray();
}

// operator

template<typename T>
Array3D<T> & Array3D<T>::operator=(const Array3D &rhs)
{
    if (this != &rhs) {

        if (n1_ != rhs.n1_ || n2_ != rhs.n2_ || n3_ != rhs.n3_) {
            DeleteArray();
            n1_ = rhs.n1_;
            n2_ = rhs.n2_;
            n3_ = rhs.n3_;
            AllocateArray();
        }

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    pv_[i][j][k] = rhs[i][j][k];
                }
            }
        }

    }

    return *this;
}

template<typename T>
Array3D<T> & Array3D<T>::operator=(const T &a)
{
    if (status_ == ArrayStatus::allocated) {

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    pv_[i][j][k] = a;
                }
            }
        }

    }

    return *this;
}

template<typename T>
inline T** Array3D<T>::operator[](const int i)
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

template<typename T>
inline const T* const * Array3D<T>::operator[](const int i) const
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// member function

template<typename T>
void Array3D<T>::Resize(const int n1, const int n2, const int n3)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        AllocateArray();
    }
}

template<typename T>
void Array3D<T>::Assign(const int n1, const int n2, const int n3, const T &a)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        AllocateArray();
    }
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                pv_[i][j][k] = a;
            }
        }
    }
}

template<typename T>
void Array3D<T>::AllocateArray()
{
    if (status_ == ArrayStatus::empty) {
        pv_       = new T**[n1_];
        pv_[0]    = new T*[n1_*n2_];
        pv_[0][0] = new T[n1_*n2_*n3_];
        for (int i = 0; i < n1_; ++i) {
            pv_[i] = pv_[0] + i*n2_;
            for (int j = 0; j < n2_; ++j) {
                pv_[i][j] = pv_[0][0] + i*n2_*n3_ + j*n3_;
            }
        }
        status_ = ArrayStatus::allocated;
    }
}

template<typename T>
void Array3D<T>::DeleteArray()
{
    if (status_ == ArrayStatus::allocated) {
        delete [] pv_[0][0];
        delete [] pv_[0];
        delete [] pv_;
        pv_ = nullptr;
        status_ = ArrayStatus::empty;
    }
}

/*
#############################################################
class Array4D
#############################################################
*/

template<typename T>
class Array4D
{
public:
    // constructor
    Array4D();
    Array4D(const int n1, const int n2, const int n3, const int n4);
    Array4D(const int n1, const int n2, const int n3, const int n4, const T &a);
    Array4D(const int n1, const int n2, const int n3, const int n4, const T *a);
    Array4D(const Array4D &rhs);

    // destructor
    ~Array4D();

    // operator
    Array4D & operator=(const Array4D &rhs);
    Array4D & operator=(const T &a);
	// typedef T value_type; 
    inline T*** operator[](const int i);
	// inline T** operator[](const int i) { return pv_[i]; };
	inline const T* const * const * operator[](const int i) const;

    // menmber function
	inline int GetDim1() const { return n1_; }
	inline int GetDim2() const { return n2_; }
    inline int GetDim3() const { return n3_; }
    inline int GetDim4() const { return n4_; }

	void Resize(const int n1, const int n2, const int n3, const int n4); 
	void Assign(const int n1, const int n2, const int n3, const int n4, const T &a);

private:
    int n1_;
    int n2_;
    int n3_;
    int n4_;
    T ****pv_;
    ArrayStatus status_;

    void AllocateArray();
    void DeleteArray();

};

// constructor

template<typename T>
Array4D<T>::Array4D() : n1_(0), n2_(0), n3_(0), n4_(0), pv_(nullptr), status_(ArrayStatus::empty) { }

template<typename T>
Array4D<T>::Array4D(const int n1, const int n2, const int n3, const int n4)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
}

template<typename T>
Array4D<T>::Array4D(const int n1, const int n2, const int n3, const int n4, const T &a)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                for (int l = 0; l < n4; ++l) {
                    pv_[i][j][k][l] = a;
                }
            }
        }
    }
}


template<typename T>
Array4D<T>::Array4D(const int n1, const int n2, const int n3, const int n4, const T *a)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                for (int l = 0; l < n4; ++l) {
                    pv_[i][j][k][l] = *a++;
                }
            }
        }
    }
}


template<typename T>
Array4D<T>::Array4D(const Array4D<T> &rhs)
    : n1_(rhs.n1_), n2_(rhs.n2_), n3_(rhs.n3_), n4_(rhs.n4_), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                for (int l = 0; l < n4_; ++l) {
                    pv_[i][j][k][l] = rhs[i][j][k][l];
                }
            }
        }
    }
}

// destructor
template<typename T>
Array4D<T>::~Array4D()
{
    DeleteArray();
}

// operator

template<typename T>
Array4D<T> & Array4D<T>::operator=(const Array4D &rhs)
{
    if (this != &rhs) {

        if (n1_ != rhs.n1_ || n2_ != rhs.n2_ || n3_ != rhs.n3_ || n4_ != rhs.n4_) {
            DeleteArray();
            n1_ = rhs.n1_;
            n2_ = rhs.n2_;
            n3_ = rhs.n3_;
            n4_ = rhs.n4_;
            AllocateArray();
        }

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    for (int l = 0; l < n4_; ++l) {
                        pv_[i][j][k][l] = rhs[i][j][k][l];
                    }
                }
            }
        }

    }

    return *this;
}

template<typename T>
Array4D<T> & Array4D<T>::operator=(const T &a)
{
    if (status_ == ArrayStatus::allocated) {

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    for (int l = 0; l < n4_; ++l) {
                        pv_[i][j][k][l] = a;
                    }
                }
            }
        }

    }

    return *this;
}

template<typename T>
inline T*** Array4D<T>::operator[](const int i)
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// template<typename T>
// inline T** Array4D<T>::operator[](const int i)
// {
// #ifdef CHECK_BOUNDS_
//     if (i < 0 || i >= n1_) {
//         throw("Array2D subscript out of bounds");
//     }
// #endif
//     return pv_[i];
// }

template<typename T>
inline const T* const * const * Array4D<T>::operator[](const int i) const
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// member function

template<typename T>
void Array4D<T>::Resize(const int n1, const int n2, const int n3, const int n4)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_ || n4 != n4_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        n4_ = n4;
        AllocateArray();
    }
}

template<typename T>
void Array4D<T>::Assign(const int n1, const int n2, const int n3, const int n4, const T &a)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_ || n4 != n4_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        n4_ = n4;
        AllocateArray();
    }
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                for (int l = 0; l < n4_; ++l) {
                    pv_[i][j][k][l] = a;
                }
            }
        }
    }
}

template<typename T>
void Array4D<T>::AllocateArray()
{
    if (status_ == ArrayStatus::empty) {
        pv_          = new T***[n1_];
        pv_[0]       = new T**[n1_*n2_];
        pv_[0][0]    = new T*[n1_*n2_*n3_];
        pv_[0][0][0] = new T[n1_*n2_*n3_*n4_];
        for (int i = 0; i < n1_; ++i) {
            pv_[i] = pv_[0] + i*n2_;
            for (int j = 0; j < n2_; ++j) {
                pv_[i][j] = pv_[0][0] + i*n2_*n3_ + j*n3_;
                for (int k = 0; k < n3_; ++k) {
                    pv_[i][j][k] = pv_[0][0][0] + i*n2_*n3_*n4_ + j*n3_*n4_ + k*n4_;
                }
            }
        }
        status_ = ArrayStatus::allocated;
    }
}

template<typename T>
void Array4D<T>::DeleteArray()
{
    if (status_ == ArrayStatus::allocated) {
        delete [] pv_[0][0][0];
        delete [] pv_[0][0];
        delete [] pv_[0];
        delete [] pv_;
        pv_ = nullptr;
        status_ = ArrayStatus::empty;
    }
}

/*
#############################################################
class Array5D
#############################################################
*/

template<typename T>
class Array5D
{
public:
    // constructor
    Array5D();
    Array5D(const int n1, const int n2, const int n3, const int n4, const int n5);
    Array5D(const int n1, const int n2, const int n3, const int n4, const int n5, const T &a);
    Array5D(const int n1, const int n2, const int n3, const int n4, const int n5, const T *a);
    Array5D(const Array5D &rhs);

    // destructor
    ~Array5D();

    // operator
    Array5D & operator=(const Array5D &rhs);
    Array5D & operator=(const T &a);
	// typedef T value_type; 
    inline T**** operator[](const int i);
    // inline T*** operator[](const int i);
	// inline T** operator[](const int i) { return pv_[i]; };
	inline const T* const * const * const * operator[](const int i) const;

    // menmber function
	inline int GetDim1() const { return n1_; }
	inline int GetDim2() const { return n2_; }
    inline int GetDim3() const { return n3_; }
    inline int GetDim4() const { return n4_; }
    inline int GetDim5() const { return n5_; }

	void Resize(const int n1, const int n2, const int n3, const int n4, const int n5); 
	void Assign(const int n1, const int n2, const int n3, const int n4, const int n5, const T &a);

private:
    int n1_;
    int n2_;
    int n3_;
    int n4_;
    int n5_;
    T *****pv_;
    ArrayStatus status_;

    void AllocateArray();
    void DeleteArray();

};

// constructor

template<typename T>
Array5D<T>::Array5D() : n1_(0), n2_(0), n3_(0), n4_(0), n5_(0), pv_(nullptr), status_(ArrayStatus::empty) { }

template<typename T>
Array5D<T>::Array5D(const int n1, const int n2, const int n3, const int n4, const int n5)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), n5_(n5), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
}

template<typename T>
Array5D<T>::Array5D(const int n1, const int n2, const int n3, const int n4, const int n5, const T &a)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), n5_(n5), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                for (int l = 0; l < n4; ++l) {
                    for (int m = 0; m < n5; ++m) {
                        pv_[i][j][k][l][m] = a;
                    }
                }
            }
        }
    }
}


template<typename T>
Array5D<T>::Array5D(const int n1, const int n2, const int n3, const int n4, const int n5, const T *a)
    : n1_(n1), n2_(n2), n3_(n3), n4_(n4), n5_(n5), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                for (int l = 0; l < n4; ++l) {
                    for (int m = 0; m < n5; ++m) {
                        pv_[i][j][k][l][m] = *a++;
                    }
                }
            }
        }
    }
}


template<typename T>
Array5D<T>::Array5D(const Array5D<T> &rhs)
    : n1_(rhs.n1_), n2_(rhs.n2_), n3_(rhs.n3_), n4_(rhs.n4_), n5_(rhs.n5_), pv_(nullptr), status_(ArrayStatus::empty)
{
    AllocateArray();
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                for (int l = 0; l < n4_; ++l) {
                    for (int m = 0; m < n5_; ++m) {
                        pv_[i][j][k][l][m] = rhs[i][j][k][l][m];
                    }
                }
            }
        }
    }
}

// destructor
template<typename T>
Array5D<T>::~Array5D()
{
    DeleteArray();
}

// operator

template<typename T>
Array5D<T> & Array5D<T>::operator=(const Array5D &rhs)
{
    if (this != &rhs) {

        if (n1_ != rhs.n1_ || n2_ != rhs.n2_ || n3_ != rhs.n3_ || n4_ != rhs.n4_ || n5_ != rhs.n5_) {
            DeleteArray();
            n1_ = rhs.n1_;
            n2_ = rhs.n2_;
            n3_ = rhs.n3_;
            n4_ = rhs.n4_;
            n5_ = rhs.n5_;
            AllocateArray();
        }

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    for (int l = 0; l < n4_; ++l) {
                        for (int m = 0; m < n5_; ++m) {
                            pv_[i][j][k][l][m] = rhs[i][j][k][l][m];
                        }
                    }
                }
            }
        }

    }

    return *this;
}


template<typename T>
Array5D<T> & Array5D<T>::operator=(const T &a)
{
    if (status_ == ArrayStatus::allocated) {

        for (int i = 0; i < n1_; ++i) {
            for (int j = 0; j < n2_; ++j) {
                for (int k = 0; k < n3_; ++k) {
                    for (int l = 0; l < n4_; ++l) {
                        for (int m = 0; m < n5_; ++m) {
                            pv_[i][j][k][l][m] = a;
                        }
                    }
                }
            }
        }

    }

    return *this;
}

template<typename T>
inline T**** Array5D<T>::operator[](const int i)
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// template<typename T>
// inline T** Array4D<T>::operator[](const int i)
// {
// #ifdef CHECK_BOUNDS_
//     if (i < 0 || i >= n1_) {
//         throw("Array2D subscript out of bounds");
//     }
// #endif
//     return pv_[i];
// }

template<typename T>
inline const T* const * const * const * Array5D<T>::operator[](const int i) const
{
#ifdef CHECK_BOUNDS_
    if (i < 0 || i >= n1_) {
        throw("Array2D subscript out of bounds");
    }
#endif
    return pv_[i];
}

// member function

template<typename T>
void Array5D<T>::Resize(const int n1, const int n2, const int n3, const int n4, const int n5)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_ || n4 != n4_ || n5 != n5_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        n4_ = n4;
        n5_ = n5;
        AllocateArray();
    }
}

template<typename T>
void Array5D<T>::Assign(const int n1, const int n2, const int n3, const int n4, const int n5, const T &a)
{
    if (n1 != n1_ || n2 != n2_ || n3 != n3_ || n4 != n4_ || n5 != n5_) {
        DeleteArray();
        n1_ = n1;
        n2_ = n2;
        n3_ = n3;
        n4_ = n4;
        n5_ = n5;
        AllocateArray();
    }
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            for (int k = 0; k < n3_; ++k) {
                for (int l = 0; l < n4_; ++l) {
                    for (int m = 0; m < n5_; ++m) {
                        pv_[i][j][k][l][m] = a;
                    }
                }
            }
        }
    }
}

template<typename T>
void Array5D<T>::AllocateArray()
{
    if (status_ == ArrayStatus::empty) {
        pv_             = new T****[n1_];
        pv_[0]          = new T***[n1_*n2_];
        pv_[0][0]       = new T**[n1_*n2_*n3_];
        pv_[0][0][0]    = new T*[n1_*n2_*n3_*n4_];
        pv_[0][0][0][0] = new T[n1_*n2_*n3_*n4_*n5_];
        for (int i = 0; i < n1_; ++i) {
            pv_[i] = pv_[0] + i*n2_;
            for (int j = 0; j < n2_; ++j) {
                pv_[i][j] = pv_[0][0] + i*n2_*n3_ + j*n3_;
                for (int k = 0; k < n3_; ++k) {
                    pv_[i][j][k] = pv_[0][0][0] + i*n2_*n3_*n4_ + j*n3_*n4_ + k*n4_;
                    for (int l = 0; l < n4_; ++l) {
                        pv_[i][j][k][l] = pv_[0][0][0][0] + i*n2_*n3_*n4_*n5_ + j*n3_*n4_*n5_ + k*n4_*n5_ + l*n5_;
                    }
                }
            }
        }
        status_ = ArrayStatus::allocated;
    }
}

template<typename T>
void Array5D<T>::DeleteArray()
{
    if (status_ == ArrayStatus::allocated) {
        delete [] pv_[0][0][0][0];
        delete [] pv_[0][0][0];
        delete [] pv_[0][0];
        delete [] pv_[0];
        delete [] pv_;
        pv_ = nullptr;
        status_ = ArrayStatus::empty;
    }
}


/*
#############################################################
using type 
#############################################################
*/

// basic type names
using Int           = int;
using Uint          = unsigned int;
using Double        = double;
using IntComplex    = std::complex<int>;
using DoubleComplex = std::complex<double>;
using Bool          = bool;
using Char          = char;
using String        = std::string;

// array type names

// I  : input (const)
// O  : output
// IO : input/output

// Array1D

using IntArray1D    = Array1D<Int>;
using IntArray1D_I  = const Array1D<Int>;
using IntArray1D_O  = Array1D<Int>;
using IntArray1D_IO = Array1D<Int>;

using UintArray1D    = Array1D<Uint>;
using UintArray1D_I  = const Array1D<Uint>;
using UintArray1D_O  = Array1D<Uint>;
using UintArray1D_IO = Array1D<Uint>;

using DoubleArray1D    = Array1D<Double>;
using DoubleArray1D_I  = const Array1D<Double>;
using DoubleArray1D_O  = Array1D<Double>;
using DoubleArray1D_IO = Array1D<Double>;

using IntComplexArray1D    = Array1D<IntComplex>;
using IntComplexArray1D_I  = const Array1D<IntComplex>;
using IntComplexArray1D_O  = Array1D<IntComplex>;
using IntComplexArray1D_IO = Array1D<IntComplex>;

using DoubleComplexArray1D    = Array1D<DoubleComplex>;
using DoubleComplexArray1D_I  = const Array1D<DoubleComplex>;
using DoubleComplexArray1D_O  = Array1D<DoubleComplex>;
using DoubleComplexArray1D_IO = Array1D<DoubleComplex>;

using BoolArray1D    = Array1D<Bool>;
using BoolArray1D_I  = const Array1D<Bool>;
using BoolArray1D_O  = Array1D<Bool>;
using BoolArray1D_IO = Array1D<Bool>;

using CharArray1D    = Array1D<Char>;
using CharArray1D_I  = const Array1D<Char>;
using CharArray1D_O  = Array1D<Char>;
using CharArray1D_IO = Array1D<Char>;

using StringArray1D    = Array1D<String>;
using StringArray1D_I  = const Array1D<String>;
using StringArray1D_O  = Array1D<String>;
using StringArray1D_IO = Array1D<String>;

// Array2D

using IntArray2D    = Array2D<Int>;
using IntArray2D_I  = const Array2D<Int>;
using IntArray2D_O  = Array2D<Int>;
using IntArray2D_IO = Array2D<Int>;

using UintArray2D    = Array2D<Uint>;
using UintArray2D_I  = const Array2D<Uint>;
using UintArray2D_O  = Array2D<Uint>;
using UintArray2D_IO = Array2D<Uint>;

using DoubleArray2D    = Array2D<Double>;
using DoubleArray2D_I  = const Array2D<Double>;
using DoubleArray2D_O  = Array2D<Double>;
using DoubleArray2D_IO = Array2D<Double>;

using IntComplexArray2D    = Array2D<IntComplex>;
using IntComplexArray2D_I  = const Array2D<IntComplex>;
using IntComplexArray2D_O  = Array2D<IntComplex>;
using IntComplexArray2D_IO = Array2D<IntComplex>;

using DoubleComplexArray2D    = Array2D<DoubleComplex>;
using DoubleComplexArray2D_I  = const Array2D<DoubleComplex>;
using DoubleComplexArray2D_O  = Array2D<DoubleComplex>;
using DoubleComplexArray2D_IO = Array2D<DoubleComplex>;

using BoolArray2D    = Array2D<Bool>;
using BoolArray2D_I  = const Array2D<Bool>;
using BoolArray2D_O  = Array2D<Bool>;
using BoolArray2D_IO = Array2D<Bool>;

// Array3D

using IntArray3D    = Array3D<Int>;
using IntArray3D_I  = const Array3D<Int>;
using IntArray3D_O  = Array3D<Int>;
using IntArray3D_IO = Array3D<Int>;

using UintArray3D    = Array3D<Uint>;
using UintArray3D_I  = const Array3D<Uint>;
using UintArray3D_O  = Array3D<Uint>;
using UintArray3D_IO = Array3D<Uint>;

using DoubleArray3D    = Array3D<Double>;
using DoubleArray3D_I  = const Array3D<Double>;
using DoubleArray3D_O  = Array3D<Double>;
using DoubleArray3D_IO = Array3D<Double>;

using IntComplexArray3D    = Array3D<IntComplex>;
using IntComplexArray3D_I  = const Array3D<IntComplex>;
using IntComplexArray3D_O  = Array3D<IntComplex>;
using IntComplexArray3D_IO = Array3D<IntComplex>;

using DoubleComplexArray3D    = Array3D<DoubleComplex>;
using DoubleComplexArray3D_I  = const Array3D<DoubleComplex>;
using DoubleComplexArray3D_O  = Array3D<DoubleComplex>;
using DoubleComplexArray3D_IO = Array3D<DoubleComplex>;

using BoolArray3D    = Array3D<Bool>;
using BoolArray3D_I  = const Array3D<Bool>;
using BoolArray3D_O  = Array3D<Bool>;
using BoolArray3D_IO = Array3D<Bool>;

// Array4D

using IntArray4D    = Array4D<Int>;
using IntArray4D_I  = const Array4D<Int>;
using IntArray4D_O  = Array4D<Int>;
using IntArray4D_IO = Array4D<Int>;

using UintArray4D    = Array4D<Uint>;
using UintArray4D_I  = const Array4D<Uint>;
using UintArray4D_O  = Array4D<Uint>;
using UintArray4D_IO = Array4D<Uint>;

using DoubleArray4D    = Array4D<Double>;
using DoubleArray4D_I  = const Array4D<Double>;
using DoubleArray4D_O  = Array4D<Double>;
using DoubleArray4D_IO = Array4D<Double>;

using IntComplexArray4D    = Array4D<IntComplex>;
using IntComplexArray4D_I  = const Array4D<IntComplex>;
using IntComplexArray4D_O  = Array4D<IntComplex>;
using IntComplexArray4D_IO = Array4D<IntComplex>;

using DoubleComplexArray4D    = Array4D<DoubleComplex>;
using DoubleComplexArray4D_I  = const Array4D<DoubleComplex>;
using DoubleComplexArray4D_O  = Array4D<DoubleComplex>;
using DoubleComplexArray4D_IO = Array4D<DoubleComplex>;

using BoolArray4D    = Array4D<Bool>;
using BoolArray4D_I  = const Array4D<Bool>;
using BoolArray4D_O  = Array4D<Bool>;
using BoolArray4D_IO = Array4D<Bool>;

// Array5D

using IntArray5D    = Array5D<Int>;
using IntArray5D_I  = const Array5D<Int>;
using IntArray5D_O  = Array5D<Int>;
using IntArray5D_IO = Array5D<Int>;

using UintArray5D    = Array5D<Uint>;
using UintArray5D_I  = const Array5D<Uint>;
using UintArray5D_O  = Array5D<Uint>;
using UintArray5D_IO = Array5D<Uint>;

using DoubleArray5D    = Array5D<Double>;
using DoubleArray5D_I  = const Array5D<Double>;
using DoubleArray5D_O  = Array5D<Double>;
using DoubleArray5D_IO = Array5D<Double>;

using IntComplexArray5D    = Array5D<IntComplex>;
using IntComplexArray5D_I  = const Array5D<IntComplex>;
using IntComplexArray5D_O  = Array5D<IntComplex>;
using IntComplexArray5D_IO = Array5D<IntComplex>;

using DoubleComplexArray5D    = Array5D<DoubleComplex>;
using DoubleComplexArray5D_I  = const Array5D<DoubleComplex>;
using DoubleComplexArray5D_O  = Array5D<DoubleComplex>;
using DoubleComplexArray5D_IO = Array5D<DoubleComplex>;

using BoolArray5D    = Array5D<Bool>;
using BoolArray5D_I  = const Array5D<Bool>;
using BoolArray5D_O  = Array5D<Bool>;
using BoolArray5D_IO = Array5D<Bool>;

}; // end namespace array

#endif /* ARRAY_HPP_ */