#ifndef VECTOR_HH
#define VECTOR_HH

template<typename T, typename U>
void Vector<T,U>::allocate() {
    rawData = new T [nx];
}

template<typename T, typename U>
Vector<T,U>::Vector(U N) : nx(N) {
    allocate();
}

template<typename T, typename U>
Vector<T,U>::Vector(const Vector<T,U> &rhs) : nx(rhs.getNx()) {
    allocate();
    for (U iX = 0; iX < nx; ++iX) {
        rawData[iX] = rhs.get(iX);
    }
}

template<typename T, typename U>
Vector<T,U>::~Vector() {
    if (rawData != 0) {
        delete[] rawData;
    }
}

template<typename T, typename U>
Vector<T,U>& Vector<T,U>::operator= (Vector<T,U> const& rhs ) {
    Vector<T,U> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, typename U>
void Vector<T,U>::swap(Vector<T,U> &rhs) {
    std::swap(rawData, rhs.rawData);
}

template<typename T, typename U>
void Vector<T,U>::print(int len) const {
    for (U iX = 0; iX < getNx(); ++iX) {
        std::cout << std::setprecision(len) << (double)get(iX) << " ";
    }
    std::cout << std::setprecision(len) <<  std::endl;
}

#endif // VECTOR_HH
