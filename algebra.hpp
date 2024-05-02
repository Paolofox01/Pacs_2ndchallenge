#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <map>
#include <array>
#include <vector>
#include <iterator>
#include <complex>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <fstream>

/*!
 * \brief Namespace containing classes for algebraic operations.
 */

namespace algebra {

    /*!
     * \brief Enumeration defining the order of compressed matrix.
     */
    enum order {
        row = 0, /*!< Row-major order */
        column = 1 /*!< Column-major order */
    };

    /*!
     * \brief Enumeration defining the type of norm to compute.
     */
    enum NormType {
        One, /*!< One-norm */
        Infinity, /*!< Infinity-norm */
        Frobenius /*!< Frobenius-norm */
    };
    
    /*!
     * \brief Functor for comparing matrix element indices based on order.
     * 
     * This struct modifies the ordering of the map between the CSR and CSC cases.
     * \tparam M Order of the matrix (row or column).
     */

     template <order M>
    struct CompOp
    { 
        /*!
         * \brief Comparison operator for matrix element indices.
         * \param l First index pair.
         * \param r Second index pair.
         * \return True if l is less than r, false otherwise.
         */
        bool operator() (std::array<std::size_t,2> const & l, std::array<std::size_t,2> const & r) const
        {
            if constexpr (M == order::row){
                // If CSR, use the default comparison operator
                return l < r;
            } else {
                // If CSC, order elements based on the column index first, then compare the row index if on the same column
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
    };

    /*!
     * \brief Class representing a sparse matrix.
     * \tparam T Type of the elements contained in the matrix.
     * \tparam N Type of compression to use (row or column).
     */
    template <typename T, order N = order::row> 
    class Matrix{
    public:
        
        /*!
         * \brief Constructor for an empty matrix.
         * \param rows Number of rows.
         * \param cols Number of columns.
         */
        Matrix (std::size_t const & rows, std::size_t const & cols) : nrows(rows), ncols(cols) { compressed = false;}
        
        /*!
         * \brief Constructor for a matrix from a .mtx file.
         * \param nomeFile Name of the .mtx file.
         */
        Matrix (const std::string &nomeFile) {

            //read the file
            std::ifstream file(nomeFile);

            if (!file.is_open()) {
                std::cerr << "couldn't open file " << nomeFile << std::endl;
            }

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

        /*!
         * \brief Resizes the matrix.
         * \param rows New number of rows.
         * \param cols New number of columns.
         */
        void resize(std::size_t const & rows, std::size_t const & cols) {

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
        }

        // operator * for the product of a matrix with a vector
        template<typename K, order M>
        friend std::vector<K> operator* (Matrix<K, M> const &, std::vector<K> const &);

        // operator * for the product of a matrix with another matrix
        template<typename K, order M>
        friend Matrix<K, M> operator* (Matrix<K, M> const &, Matrix<K, M> const &);

        // & operator () to modify the elements of the matrix
        /*!
        * \brief Accesses and modifies elements of the matrix.
        * \param i Row index.
        * \param j Column index.
        * \return Reference to the matrix element at the specified position.
        */
        auto &
        operator()(std::size_t const & i, std::size_t const & j) {

            //control if the indexes are inside the limits of the matrix
            
            if(i < nrows && j < ncols ) {

                //if compressed return the address to the element if !=0, else print an error( can't modify a 0 element if compressed)          
                
                if (compressed) {
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
        /*!
 * \brief Accesses elements of the matrix.
 * \param i Row index.
 * \param j Column index.
 * \return The value of the matrix element at the specified position, or 0 if not present but within the limits.
 */
        auto
        operator()(std::size_t const & i, std::size_t const & j) const
        {
            //very similar to the & operator () version, but this time return 0 if the element is not present but inside the limits
            if(i < nrows && j < ncols ) {
                
                if (compressed) {
                    
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
            
        //this function trasforms the matrix from uncopressed to compressed
        /*!
 * \brief Compresses the matrix representation.
 */
        void compress() {

            if(compressed) {
                std::cerr << "already compressed" << std::endl;
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

        //this function trasforms the matrix from uncopressed to compressed

        /*!
 * \brief Decompresses the matrix representation.
 */
        void uncompress() {
            
            if(!compressed) {
                std::cerr << "already uncompressed" << std::endl;
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
            inner_i.clear();
            outer_i.clear();
            values.clear();

            //update the compressed status
            compressed = 0;

            return;
        }

        //returns whether the matrix is in compressed form 
        /*!
 * \brief Checks if the matrix is in compressed form.
 * \return True if the matrix is compressed, false otherwise.
 */
        bool is_compressed() {
            return compressed;
        }


        //returns the number of rows in the matrix
                /*!
        * \brief Getter for nrows
        */
        auto rows() const {return nrows };

        //returns the number of columns in the matrix
                /*!
        * \brief Getter for ncols
        */
        auto cols() const {return ncols };


        // function that computes the norm of a matrix depending on the template 
        /*!
        * \brief Computes the norm of the matrix.
        * \tparam tipo Type of norm to compute (One, Infinity, or Frobenius).
        * \return The computed norm value.
        */
        template <NormType tipo = NormType::One>
        double norm() {
            
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
                        totali[outer_i[i]] += values[i]; 
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
                        totali[outer_i[i]] += values[i]; 
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



        private:

        //map that contains an array with the coordinates and the corrisponding value, ordered based on CompOp
        std::map<std::array<std::size_t,2>,T, CompOp<N> > numbers;
        
        //vector that contains the inner indexes of the matrix saved in compressed form
        std::vector<std::size_t> inner_i;
        
        //vector that contains the outer indexes of the matrix saved in compressed form
        std::vector<std::size_t> outer_i;
        
        //vector that contains the values inside the matrix saved in compressed form
        std::vector<T> values;
        
        //number of rows
        std::size_t nrows;
        
        //number of columns
        std::size_t ncols;
        
        //true if in compressed form, false otherwise
        bool compressed;
    };




    //operator to compute the product of a matrix and a vector
    /*!
 * \brief Computes the product of a matrix with a vector.
 * \param Mat Matrix operand.
 * \param vec Vector operand.
 * \return The resulting vector.
 */
    template <typename T, order N>
    std::vector<T> operator* (Matrix<T, N> const & Mat, std::vector<T> const & vec) {
        
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


    //operator to compute the product of two matrixes
    /*!
 * \brief Computes the product of two matrices.
 * \param Mat_l Left matrix operand.
 * \param Mat_r Right matrix operand.
 * \return The resulting matrix.
 */
    template <typename T, order N>
    Matrix<T, N> operator* (Matrix<T, N> const & Mat_l, Matrix<T, N> const & Mat_r) {
        

        if( Mat_l.ncols != Mat_r.nrows) {
            std::cerr << "wrong dimensions" << std::endl;
            return algebra::Matrix<T, N>(0, 0);
        }    

        if( Mat_l.compressed + Mat_r.compressed == 1) {
            std::cerr << "not in the same form" << std::endl;
            return algebra::Matrix<T, N>(0, 0);
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