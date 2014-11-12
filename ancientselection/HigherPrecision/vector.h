#ifndef VECTOR_H
#define VECTOR_H


template<typename T, typename U>
class Vector {
public :
    Vector() { };
    Vector(U N);
    Vector(const Vector<T,U> &rhs);
    ~Vector();
    
    Vector<T,U>& operator=(Vector<T,U> const& rhs);
    void swap(Vector<T,U> &rhs);
public :    
    U getNx() const {
        return nx;
    }
    
    T *get() {
        return rawData;
    }
    
    T & get(U iX) {
        assert(iX < nx);
        return rawData[iX];
    }
    
    T const& get(U iX) const {
        assert(iX < nx);
        return rawData[iX];
    }
    void print(int len = 15) const;
private :
    void allocate();
    
private :
    U nx;
    T *rawData;
};

#endif // VECTOR_H
