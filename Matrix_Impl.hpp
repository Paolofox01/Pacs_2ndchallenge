#ifndef MATRIXIMPL_H
#define MATRIXIMPL_H
// clang-format off
#include "algebra.hpp"

namespace algebra {

    template <order M> 
    bool CompOp<M>::operator() (std::array<std::size_t,2> const & l, std::array<std::size_t,2> const & r) const 
        {
            if constexpr (M == order::row){
                // If CSR, use the default comparison operator
                return l < r;
            } else {
                // If CSC, order elements based on the column index first, then compare the row index if on the same column
                //@note A trick using std::tie can be used to simplify the code
                // return std::tie(l.at(1), l.at(0)) < std::tie(r.at(1), r.at(0));
                if (l.at(1) < r.at(1))
                    return true;
                else if (l.at(1) > r.at(1))
                    return false;
                else {
                    if (l.at(0) < r.at(0))
                        return true;
                    else return false;
                }
            }
        }

    // Constructor for a matrix from a .mtx file
    template<typename T, order N>
    Matrix<T, N>::Matrix(const std::string &nomeFile)  {

            //read the file
            std::ifstream file(nomeFile);

            if (!file.is_open()) {
                std::cerr << "couldn't open file " << nomeFile << std::endl;
            }
 //@note It is better to read a line with getline and then process the string also to eliminate comment lines

            while (file.peek() == '%') file.ignore(2048, '\n');

            std::string riga;

            std::getline(file, riga);

            std::istringstream ss(riga);
            //the first row contains the dimensions of the matrix
            ss >> nrows >> ncols;

            std::size_t rigaIndice, colonnaIndice;
            double valore;

            //fill the map with the values
            while (file >> rigaIndice >> colonnaIndice >> valore) {
                numbers[{rigaIndice - 1, colonnaIndice -1}] = valore; 
            }

            compressed = false;

        }

    // Resizes the matrix
    template<typename T, order N>
    void Matrix<T, N>::resize(std::size_t const &rows, std::size_t const &cols)  {

            // resize works with uncompressed matrixes

            if(compressed){
            uncompress();}

            // if the resized matrix is smaller than the original matrix erase the elements outside the limits 

            if(rows < nrows || cols < ncols){
                
                for (auto it = numbers.begin(); it != numbers.end(); )
                {
                    if (it->first.at(0) >= rows || it->first.at(1) >= cols) {
                        
                        it = numbers.erase(it);

                    } else {
                        it++;
                    }
                }
                
            }
            
            // update nrows and ncols

            nrows = rows;
            ncols = cols;
            //@note OK, but if the matrix was originally compresses I would recompress it again. To avoid confusion from the user.
        }

    // & operator () to modify the elements of the matrix
    template<typename T, order N>
    auto & Matrix<T, N>::operator()(std::size_t const & i, std::size_t const & j) {

            //control if the indexes are inside the limits of the matrix
            
            if(i < nrows && j < ncols ) {

                //if compressed return the address to the element if !=0, else print an error( can't modify a 0 element if compressed)          
                
                if (compressed) {
                    //@note I don't understand what you are doing. If the matrix is rowwise ordered
                    // you should look for for each row i all the elements in inner_i[outer_i[i]] to inner_i[outer_i[i+1\] 
                    // and check if the column index is equal to j. If the matrix is column oriented the role of inner_i and outer_i is reversed.
                    for (auto k = inner_i[i*(1-N) + j * N]; k < inner_i[i*(1-N) + j * N + 1]; k++){
                        if(outer_i[k] == i * N + j * (1-N))
                            return values[k]; 
                    }
                    std::cerr << "attempt to modify a zero element of a compressed matrix" << std::endl;                 
                }
                    
                //if uncompressed just return the address of the element in the map corresponding to the to the coordinates
                
                else {
                    return numbers[{i, j}];
                }
            }
            else
            {
                std::cerr << "out of bounds" << std::endl;
            }

        std::exit(EXIT_FAILURE);
        }


    //const oprator () to read the elements from the matrix
    template<typename T, order N>
    auto Matrix<T, N>::operator()(std::size_t const & i, std::size_t const & j) const
        {
            //very similar to the & operator () version, but this time return 0 if the element is not present but inside the limits
            if(i < nrows && j < ncols ) {
                
                if (compressed) {
                    //@note same comment as before
                    for (auto k = inner_i[i*(1-N) + j * N]; k < inner_i[i*(1-N) + j * N + 1]; k++)
                        if(outer_i[k] == i*N + j * (1-N))
                            return values[k]; }
                else{
                    if(numbers.find({i, j}) != numbers.end())
                        return numbers.at({i, j});
                }
            }
            else
            {
                std::cerr << "out of bounds" << std::endl;
            }
            
            return 0;

        }

    // Compresses the matrix representation
    template<typename T, order N>
    void Matrix<T, N>::compress()  {

            if(compressed) {
                std::cerr << "already compressed" << std::endl;//@note it is not necessarily an error
                return;
            }
            
            // allocate the right mmempory for the vectors

            outer_i.resize(numbers.size());
            
            values.resize(numbers.size());
            
            int k = 0;
            
            auto iniz = numbers.begin();

            //insert the values in the vectors using the right order type
            
            if constexpr (N == order::row) {
                inner_i.resize(nrows + 1);
                    
                for (auto it = iniz; it != numbers.end(); it++) {
                    outer_i[k] = it->first.at(1);
                    values[k] = it->second;
                    k++;
                }
                
                for(std::size_t i = 0; i < nrows+1; i++) {
//@note no need to use distance. You can exploit the fact theat the map is ordered, so you can extract row by row (or col by col if the matrix is column oriented) the elements
                    inner_i[i] = std::distance(iniz, numbers.lower_bound(std::array<std::size_t,2>{i, 0}));         
                    
            } 
            } else {
                inner_i.resize(ncols + 1);
                
                for (auto it = iniz; it != numbers.end(); it++) {
                    outer_i[k] = it->first.at(0);
                    values[k] = it->second;
                    k++;
                }
                
                for(std::size_t i = 0; i < ncols+1; i++) {

                    inner_i[i] = std::distance(iniz, numbers.lower_bound(std::array<std::size_t,2>{0, i}));         
                    
            }

            }

            // free the memory that contains the map

            numbers.clear();

            //update the compressed status

            compressed = 1;

            return;
        }

    // Decompresses the matrix representation
    template<typename T, order N>
    void Matrix<T, N>::uncompress()  {
            
            if(!compressed) {
                std::cerr << "already uncompressed" << std::endl;//@note not necessarily an error
                return;
            }

            // fill the map

            std::size_t k = 0;
            for (size_t i = 0; i < values.size(); i++)
            {
                while( inner_i[k+1] <= i )
                k++;

                if constexpr (N == order::row) { 
                    numbers[{k, outer_i[i]}] = values[i];
                } else {
                    numbers[{outer_i[i], k}] = values[i];
                }

            }

            //free the memory
            //@note you are not freeing the memory with clear() on a vector. You should use shrink_to_fit() after the clear() to free the memory
            inner_i.clear();
            outer_i.clear();
            values.clear();

            //update the compressed status
            compressed = 0;

            return;
        }

    // Computes the norm of the matrix
    template<typename T, order N>
    template <NormType tipo>
    double Matrix<T, N>::norm(){ //@note This method should be const. It does not change the state of the matrix!
            //@note Nice having this method working also for the uncompressed case. Yet normally for this high lever
            // method I would expect to have a compressed matrix. 
            //apply the formulas to compute the norm to the matrix based on the type of norm and matrix ordering
            if constexpr (N == order::row) {
            if constexpr (tipo == NormType::One) {
                 std::vector<double> totali(ncols, 0.0);
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            totali[key.at(1)] += std::abs(value);
                        }
                    return *std::max_element(totali.cbegin(), totali.cend());
                    
                } else {
                    for (std::size_t i=0; i < outer_i.size(); i++) {
                        totali[outer_i[i]] += std::abs(values[i]); 
                    }

                    return *std::max_element(totali.cbegin(), totali.cend());
                }

            } else if constexpr (tipo == NormType::Infinity) {
                 std::vector<double> totali(nrows, 0.0);
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            totali[key.at(0)] += std::abs(value);
                        }
                    return *std::max_element(totali.cbegin(), totali.cend());
                    
                } else {
                    for (std::size_t i=0; i < inner_i.size(); i++) {
                        for(std::size_t j = inner_i[i]; j < inner_i[i+1]; j++)
                        totali[i] +=  std::abs(values[j]); 
                    }

                    return *std::max_element(totali.cbegin(), totali.cend());
                }


            } else {
                double res = 0.0;
                
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            //@note don't use pow just to elevete to 2. Use value * value (less expensive)
                            res += std::pow(std::abs(value), 2);
                        }
                    return std::sqrt(res);
                    
                } else {
                    for (auto i : values) {
                        res += std::pow(std::abs(i), 2);//@note same comment as before. Moreover why taking abs when you are squaring?
                    }

                    return std::sqrt(res);
                } 
                } 
            } else {
                if constexpr (tipo == NormType::One) {
                 std::vector<double> totali(ncols, 0.0);
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            totali[key.at(1)] += std::abs(value);
                        }
                    return *std::max_element(totali.cbegin(), totali.cend());
                    
                } else {
                    for (std::size_t i=0; i < inner_i.size(); i++) {
                        for(std::size_t j = inner_i[i]; j < inner_i[i+1]; j++)
                        totali[i] +=  std::abs(values[j]);
                    }

                    return *std::max_element(totali.cbegin(), totali.cend());
                }

            } else if constexpr (tipo == NormType::Infinity) {
                 std::vector<double> totali(nrows, 0.0);
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            totali[key.at(0)] += std::abs(value);
                        }
                    return *std::max_element(totali.cbegin(), totali.cend());
                    
                } else {

                     for (std::size_t i=0; i < outer_i.size(); i++) {
                        totali[outer_i[i]] += std::abs(values[i]); 
                    }

                    return *std::max_element(totali.cbegin(), totali.cend());
                }


            } else {
                double res = 0.0;
                
                if(!compressed) {
                    for (auto [ key , value ] : numbers)
                        {
                            res += std::pow(std::abs(value), 2);
                        }
                    return std::sqrt(res);
                    
                } else {
                    for (auto i : values) {
                        res += std::pow(std::abs(i), 2);
                    }

                    return std::sqrt(res);
                } 
                }
            }
        }

    // Computes the product of a matrix with a vector
    template<typename T, order N>
    std::vector<T> operator*(Matrix<T, N> const &Mat, std::vector<T> const &vec) {
        
        //check the dimensions
        if( vec.size() != Mat.ncols) {
            std::cerr << "wrong dimensions" << std::endl;
            return std::vector<T>{};
        }    
        
        //compute the product based on how the matrix is saved
        if(Mat.compressed) {
            if constexpr (N == order::row) {
            std::vector<T> res(Mat.nrows, 0);

            for (size_t i = 0; i < Mat.nrows; i++)
            {
                
                for (size_t j = Mat.inner_i[i]; j < Mat.inner_i[i+1]; j++)
                {
                    res[i] += Mat.values[j] * vec[Mat.outer_i[j]];
                }

            }            

            return res;
        } else {
                std::vector<T> res(Mat.nrows, 0.0);

                for (size_t i = 0; i < Mat.ncols; i++)
                {
                    
                    for (size_t j = Mat.inner_i[i]; j < Mat.inner_i[i+1]; j++)
                    {
                        res[Mat.outer_i[j]] += Mat.values[j] * vec[i];
                    }
                    
                }


                return res;                

            }

        } else { 

            std::vector<T> res(Mat.nrows, 0);

                for (auto it = Mat.numbers.begin(); it != Mat.numbers.end(); it++){
                    
                    res[it->first.at(0)] += it->second * vec[it->first.at(1)];                

                }
            
            return res;

        }

    return std::vector<T>{};

    }


    // Computes the product of two matrices
    template<typename T, order N>
    Matrix<T, N> operator*(Matrix<T, N> const &Mat_l, Matrix<T, N> const &Mat_r)  {
        

        if( Mat_l.ncols != Mat_r.nrows) {
            std::cerr << "wrong dimensions" << std::endl;
            return algebra::Matrix<T, N>(0, 0);
        }    

        //to compute the multiplication with one matrix in compressed form and the other in dynamic form we use them both in compressed form

        if( Mat_l.compressed + Mat_r.compressed == 1) {
            if (Mat_l.compressed)
            { 
                Matrix<T,N> Mat_new = Mat_l;
                Mat_new.uncompress();
                return Mat_new * Mat_r;
            } else {
                Matrix<T,N> Mat_new = Mat_r;
                Mat_new.uncompress();
                return Mat_l * Mat_new;
            }
            
        }    
        
        if(Mat_l.compressed && Mat_r.compressed ) {
            if constexpr (N == order::row) {

                algebra::Matrix<T, N> res(Mat_l.nrows, Mat_r.ncols);

                for (size_t i = 0; i < Mat_l.nrows; i++)
                {
                    for (size_t j = Mat_l.inner_i[i]; j < Mat_l.inner_i[i+1]; j++)
                    {
                            for (size_t k = Mat_r.inner_i[Mat_l.outer_i[j]]; k < Mat_r.inner_i[Mat_l.outer_i[j] + 1]; k++) {
                            
                                res( i, Mat_r.outer_i[k]) += Mat_l.values[j] * Mat_r.values[k];

                            }
                        
                    }
                }
                return res;

            } else {

                algebra::Matrix<T, N> res(Mat_l.nrows, Mat_r.ncols);

                for (size_t i = 0; i < Mat_r.ncols; i++)
                {
                    for (size_t j = Mat_r.inner_i[i]; j < Mat_r.inner_i[i+1]; j++)
                    {
                            for (size_t k = Mat_l.inner_i[Mat_r.outer_i[j]]; k < Mat_l.inner_i[Mat_r.outer_i[j] + 1]; k++) {
                            
                                res(Mat_l.outer_i[k],  i) += Mat_r.values[j] * Mat_l.values[k];

                            }
                        
                    }
                }
                return res;
                

            }

        } else { 

            algebra::Matrix<T, N> res(Mat_l.nrows, Mat_r.ncols);

            for (auto [ key1 , value1 ] : Mat_l.numbers) {
                for (auto [ key2 , value2 ] : Mat_r.numbers) {
                    if (key1.at(1) == key2.at(0))
                    {
                        res(key1.at(0), key2.at(1)) += value1 * value2;
                    }
                    

            }

            }
            
            return res;

        }

    return algebra::Matrix<T, N>(0, 0);

    }
}


#endif

