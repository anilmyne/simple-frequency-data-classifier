# include <iostream>
# include <ostream>
# include <string>
# include <cstring>
# include <vector>
# include <numeric>
# include <algorithm>
# include <cmath>

#ifndef ALGORITHMS_HPP_INCLUDED
#define ALGORITHMS_HPP_INCLUDED

typedef std::vector<int> VECTOR_INT;
typedef std::vector< std::vector <int> > VECTOR_INT_2D;
typedef std::vector<long> VECTOR_LONG;
typedef std::vector< std::vector <long> > VECTOR_LONG_2D;
typedef std::vector<float> VECTOR_FLOAT;
typedef std::vector< std::vector <float> > VECTOR_FLOAT_2D;
typedef std::vector<double> VECTOR_DOUBLE;
typedef std::vector< std::vector <double> > VECTOR_DOUBLE_2D;

// atl NAMESPACE
namespace atl{
	const double PI = 2 * std::acos(0.0);
	const double V_GRAD = 1.633123935319537e+16; // vertical gradient

	float nround(double d, int r=3);
	float fround(double d, int r=3);
	float cround(double d, int r=3);

	bool isnum(std::string data);
	bool is_lower(std::string s);
	bool is_upper(std::string s);
	std::string to_upper(std::string data);
	std::string to_lower(std::string data);
	void r_to_upper(std::string *data);
	void r_to_lower(std::string *data);

	template<typename T>
	long long sum(std::vector<T> v);

	class point;
    std::vector<point> vec2points(std::vector<std::vector<int> > intpoints);
	template <class T>
	bool isequal(std::vector<T> v1, std::vector<T> v2);

	template <typename T>
	T max(std::vector<T> v){
		T val =  *std::max_element(v.begin(), v.end());
		return val;
	}
	template <typename T>
	T min(std::vector<T> v){
		T val =  *std::min_element(v.begin(), v.end());
		return val;
	}
	template <typename T>
	int indexOfmax(std::vector<T> v){
		return std::max_element(v.begin(),v.end()) - v.begin();
	}
	template <typename T>
	int indexOfmin(std::vector<T> v){
		return std::min_element(v.begin(),v.end()) - v.begin();
	}
	template <typename T>
	std::vector<int> minmax(std::vector<T> v){
		auto minmax = std::minmax_element(v.begin(), v.end());
		return minmax;
	}

	template <class T>
	int indexOf(T item, std::vector<T>);
	template <class T> // also long cast
	std::vector<int> indicesOf(T item, std::vector<T>);
	std::string readfile(std::string filename);
	std::string writefile(std::string filename, std::string data);
	template <class T>
	bool isequalarray(std::vector<T> v1, std::vector<T> v2);

	template <class T, typename U>
	std::vector<T> slice(std::vector<T> v, U start, U end=0);
	template <class T, typename U>
	void r_slice(std::vector<T> *v, U start, U end=0);
	template <class T>
	void addvector(std::vector<T> *v1, std::vector<T> v2);
	template <class T>
	void printvector(std::vector<T> v, int breaksize=0){
	    std::cout << "{";
	    for (int i = 0, counter=0; i < v.size(); i++){
            std::cout << v.at(i);
            if (i < v.size() - 1){ std::cout << ", "; }
            if (breaksize != 0){
                counter++;
                if (counter == breaksize){
                    std::cout << "\n";
                    counter = 0;
                }
            }
        }
	    std::cout << "}";
	};
	template <class T>
	void printvector(std::vector< std::vector<T> > v, int breaksize=0){
	    std::cout << "{";
	    for (int i = 0, counter=0; i < v.size(); i++){
            printvector(v[i]);
            if (i < v.size() - 1){ std::cout << ", "; }
            if (breaksize != 0){
                counter++;
                if (counter == breaksize){
                    std::cout << "\n";
                    counter = 0;
                }
            }
        }
	    std::cout << "}";
	};
	template<class T>
	int compare(T a, T b);
	template<class T>
	int compare(T *a, T *b, size_t la, size_t lb);
	template<class T>
	int compare(std::vector<T> a, std::vector<T> b);

	namespace priv{ // private anonymous sector
        template <class T>
        size_t qpartition(std::vector<T> &items, size_t left, size_t right, size_t pivot_index){
            T pcache = items[pivot_index], pivot = items[pivot_index];
            items.erase(items.begin() + pivot_index);
            items.insert(items.begin()+right, pcache);
            size_t store_index = left;
            for (size_t i=left; i < right; i++){ // range is non-inclusive of ending value
                if (atl::compare(items[i], pivot) == -1){
                    pcache = items[i];
                     items.erase(items.begin() + i);
                    items.insert(items.begin() + store_index, pcache);
                    store_index++;
                }
            }
            pcache = items[right];
             items.erase(items.begin() + right);
            items.insert(items.begin() + store_index, pcache); // move pivot to correct position
            return store_index;
        }
    }

    template <class T>
    void qsort(std::vector<T> &items, size_t left, size_t right){
        if (left < right){
            size_t pivot_index = left;
            size_t pivot_new_index = priv::qpartition(items, left, right, pivot_index);
            qsort(items, left, pivot_new_index - 1);
            qsort(items, pivot_new_index + 1, right);
        }
    }
    double crossproduct(point O, point A, point B);
    std::vector<atl::point> convex_hull(std::vector<atl::point> P);
    template <class T>
    int binarySearch(T arr[], T item, int p, int r);
    template <class T>
    int binarySearch(std::vector<T> arr, T item);
    template <class T>
    int upperindexofappend(std::vector<T> v, T item);
}
// END OF atl

typedef std::vector<atl::point> VECTOR_POINT;
typedef std::vector< std::vector <atl::point> > VECTOR_POINT_2D;


template<typename T>
long long atl::sum(std::vector<T> v){
	return std::accumulate(v.begin(), v.end(), 0);
}

class atl::point{
	public:
		point(int x, int y){
			x_ = x; y_ = y;
		}
		point(){
			x_ = 0; y_ = 0;
		}
		int x() const {return x_;}
		int y() const { return y_;}
		void x(int x){x_ = x;}
		void y(int y){y_ = y;}
		void show(){std::cout << "(" << x_ << ", " << y_ << ")";}
		void setprintprec(int p){printprec = p;};
		std::string toString();
		bool operator ==(const point &pt) const;
		bool operator !=(const point &pt) const;
		bool operator < (const point &pt) const;
		bool operator > (const point &pt) const;
		int operator [](int index) const;
		//friend std::ostream& operator <<(std::ostream& os, const point& pt);
		template <typename Rt>
		Rt distance(point pt);
	private:
		int printprec = 4;
		int x_, y_;


};

bool atl::point::operator ==(const atl::point& pt) const{
	if((x() == pt.x()) and (y() == pt.y())){
		return true;
	}return false;
}

bool atl::point::operator !=(const atl::point& pt) const{
	if((x() != pt.x()) or (y() != pt.y())){
		return true;
	}return false;
}

int atl::point::operator [](int index) const{
	if (index == 0){ return this->x(); }
	else if (index == 1) {return this->y(); }
	else {std::cout << "\npoint has index in range(0, 2)"; return -40000;}
}

std::ostream& operator << (std::ostream& os, const atl::point& pt){
    std::cout << "(" << pt.x() << ", " << pt.y() << ")";
    return os;
}

bool atl::point::operator < (const point &pt) const{
    if (this->x() == pt.x()){
        if (this->y() < pt.y()){ return true; }
    }else{
        if (this->x() < pt.x()){ return true; }
    }
    return false;
}

bool atl::point::operator > (const point &pt) const{
    if (this->x() == pt.x()){
        if (this->y() > pt.y()){ return true; }
    }else{
        if (this->x() > pt.x()){ return true; }
    }
    return false;
}

template <class T>
std::vector<int> atl::indicesOf(T item, std::vector<T> array){
    // returns all indices of item in array, array length, item in array//
    std::vector<int> indices;
    if (array.size() == 0 or indexOf(item, array) == -1)
        return indices;
    int count = std::count(array.begin(), array.end(), item);
    for (int i = 0; i < array.size(); i++){
        if (array[i] == item){ indices.push_back(i); }
        if (indices.size() == count){ break; }
	}
    return indices;
}

template <class T, typename U>
std::vector<T> atl::slice(std::vector<T> v, U start, U end){
	std::vector<T> tmp;
	if (end == 0){
		end = v.size();
	}else if (end < 0){
		end += v.size();
	}
	U i = start;
	while (i < end){
		tmp.push_back(v[i]);
		i++;
	}
	return tmp;
}

template <class T, typename U>
void atl::r_slice(std::vector<T> *v, U start, U end){
	std::vector<T> tmp;
	if (end == 0){
		end = v->size();
	}else if (end < 0){
		end += v->size();
	}
	U i = start;
	while (i < end){
		tmp.push_back(v->at(i));
		i++;
	}
	v->clear();
	for (int i = 0; i < tmp.size(); i++)
		v->push_back(tmp.at(i));
}

float atl::nround(double d, int r){ // math.h
	return roundf(d * (pow(10, r))) / pow(10, r);
}

float atl::fround(double d, int r){ // round the rounded down or trunc
	return floorf(d * (pow(10, r))) / pow(10, r);
}

float atl::cround(double d, int r){ // round the rounded up
	return ceilf(d * (pow(10, r))) / pow(10, r);
}

bool atl::isnum(std::string data){
	for(int i = 0; i < (int) data.size(); i++){
		char c = (char) data.at(i);
		if (not std::isdigit(c)){
			return false;
		}
	}
	return true;
}

std::string atl::to_upper(std::string data){
	for(int i = 0; i < data.size(); i++){
		data.at(i) = std::toupper(data.at(i));
	}
	return data;
}

std::string atl::to_lower(std::string data){
	for(int i = 0; i < data.size(); i++){
		data.at(i) = std::tolower(data.at(i));
	}
	return data;
}

bool atl::is_lower(std::string s){
	for (int i = 0; i < s.size(); i++){
		if ( not islower(s.at(i)))
			return false;
	}
	return true;
}

template <class T>
void addvector(std::vector<T> *v1, std::vector<T> v2){
    v1->insert(v1->end(), v2->begin(), v2->end());
}

template<class T>
int atl::compare(T a, T b){
    if (a == b){ return 0; }
    else if (a < b){ return -1; }
    else{ return 1; }
}

template<class T>
int atl::compare(T *a, T *b, size_t la, size_t lb){
    if (la != lb){
        if (la > lb){ return 1; }
        else{ return -1; }
    }else{
        for (size_t i=0; i < la; i++){
            if (atl::compare<T>(*(a + i), *(b + i)) != 0){
                return atl::compare<T>(*(a + i), *(b + i));
            }
        }
    }
}
template<class T>
int atl::compare(std::vector<T> a, std::vector<T> b){
    if (a.size() != b.size()){
        if (a.size() > b.size()){ return 1; }
        else{ return -1;}
    }else{
        for (size_t i=0; i < a.size(); i++){
            if (atl::compare<T>(a[i], b[i]) != 0){
                return atl::compare<T>(a[i], b[i]);
            }
        }
    }
    return -1; //
}

std::vector<atl::point> atl::vec2points(std::vector<std::vector<int> > intpoints){
    std::vector<atl::point> ret;
    for (int i=0; i < intpoints.size(); i++){
        ret.push_back(atl::point(intpoints[i][0], intpoints[i][1]));
    }
    return ret;
}

double atl::crossproduct(atl::point O, atl::point A, atl::point B){
    // wikipedia :Convex hull/Monotone chain
	return double((A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]));
}

// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
std::vector<atl::point> atl::convex_hull(std::vector<atl::point> P){
    // wikipedia :Convex hull/Monotone chain
	size_t n = P.size(), k = 0;
	if (n <= 3) return P;
	std::vector<atl::point> H(2*n);
	sort(P.begin(), P.end());
	// Build lower hull
	for (size_t i = 0; i < n; ++i) {
		while (k >= 2 && atl::crossproduct(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
	// Build upper hull
	for (size_t i = n-1, t = k+1; i > 0; --i) {
		while (k >= t && atl::crossproduct(H[k-2], H[k-1], P[i-1]) <= 0) k--;
		H[k++] = P[i-1];
	}

	H.resize(k-1);
	return H;
}

template <class T>
int atl::binarySearch(T arr[], T item, int p, int r){
    int mid;
    if (r == 0) { return -1; }
    while(p <= r){
        mid = (p + r) / 2;
        if (mid == arr.size()){ break; }
        if (arr[mid] == item){ return mid; }
        if (arr[mid] > item){ r = mid -1; }
        else{ p= mid + 1; }
    }
    return -1;
}

template <class T>  // use compare
int atl::binarySearch(std::vector<T> arr, T item){
    int p = 0, r = arr.size(), mid;
    if (r == 0) { return -1; }
    while(p <= r){
        mid = (p + r) / 2;
        if (mid == arr.size()){ break; }
        if (arr[mid] == item){ return mid; }
        if (arr[mid] > item){ r = mid -1; }
        else{ p = mid + 1; }
    }
    return -1;
}

template <class T>
int atl::upperindexofappend(std::vector<T> v, T item){
    // index where to append inorder to maintain sorted (log2N)
    int p = 0, r = v.size(), mid;
    if (r == 0) { return 0; }
    while(p <= r){
        mid = (p + r) / 2;
        if (mid == v.size()){ break; }
        if (v[mid] > item){ r = mid - 1; }
        else{ p = mid + 1; }
    }
    while(1){
        if (mid < v.size()){
            if(v[mid] < item || v[mid] == item){ mid++; }
            else{ break;}
        }else{
            break;
        }
    }
    return mid;
}

#endif
