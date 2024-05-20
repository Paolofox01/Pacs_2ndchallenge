#ifndef ALGEBRA_H
#define ALGEBRA_H
// clang-format off
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
        bool operator() (std::array<std::size_t,2> const & l, std::array<std::size_t,2> const & r) const;
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
        Matrix (std::size_t const & rows, std::size_t const & cols) : nrows(rows), ncols(cols) { compressed = false;}//@note why not initialise compress in the initialization list??
        
        //@note not needed. You can use the syntethised one. Don't define copy or assignement operator if the syntethised ones are ok
        Matrix (Matrix<T, N> const & Mat) : numbers(Mat.numbers), inner_i(Mat.inner_i), outer_i(Mat.outer_i), values(Mat.values), nrows(Mat.nrows), ncols(Mat.ncols), compressed(Mat.compressed) {}

        /*!
         * \brief Constructor for a matrix from a .mtx file.
         * \param nomeFile Name of the .mtx file.
         */
        Matrix (const std::string &nomeFile);

        /*!
         * \brief Resizes the matrix.
         * \param rows New number of rows.
         * \param cols New number of columns.
         */
        void resize(std::size_t const & rows, std::size_t const & cols);

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
        operator()(std::size_t const & i, std::size_t const & j);
        
        //const oprator () to read the elements from the matrix
        /*!
 * \brief Accesses elements of the matrix.
 * \param i Row index.
 * \param j Column index.
 * \return The value of the matrix element at the specified position, or 0 if not present but within the limits.
 */
        auto
        operator()(std::size_t const & i, std::size_t const & j) const;
            
        //this function trasforms the matrix from uncopressed to compressed
        /*!
 * \brief Compresses the matrix representation.
 */
        void compress();

        //this function trasforms the matrix from uncopressed to compressed

        /*!
 * \brief Decompresses the matrix representation.
 */
        void uncompress();

        //returns whether the matrix is in compressed form 
        /*!
 * \brief Checks if the matrix is in compressed form.
 * \return True if the matrix is compressed, false otherwise.
 */
        bool is_compressed() { return compressed; }


        //returns the number of rows in the matrix
                /*!
        * \brief Getter for nrows
        */
        auto rows() const {return nrows; };

        //returns the number of columns in the matrix
                /*!
        * \brief Getter for ncols
        */
        auto cols() const {return ncols; };


        // function that computes the norm of a matrix depending on the template 
        /*!
        * \brief Computes the norm of the matrix.
        * \tparam tipo Type of norm to compute (One, Infinity, or Frobenius).
        * \return The computed norm value.
        */
        template <NormType tipo = NormType::One>
        double norm();

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
    std::vector<T> operator* (Matrix<T, N> const & Mat, std::vector<T> const & vec);

    //operator to compute the product of two matrixes
    /*!
 * \brief Computes the product of two matrices.
 * \param Mat_l Left matrix operand.
 * \param Mat_r Right matrix operand.
 * \return The resulting matrix(if the input matrixes are both in compressed form the result will be in compressed form too, else it will be in dynamic form).
 */
    template <typename T, order N>
    Matrix<T, N> operator* (Matrix<T, N> const & Mat_l, Matrix<T, N> const & Mat_r);


}

#endif