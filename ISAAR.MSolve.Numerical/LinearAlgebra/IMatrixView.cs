using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    /*
     * 1) Determinants, inverses, factorizations, eigenvalues, SVDs, etc have complexity T(n)>O(n^2) and should not be 
     * put here. 2x2 and 3x3 matrices have analytic formulas for these operations and should expose them on their 
     * concrete classes. Not sure about matrix-matrix multiplication which is also O(n^3) (without divide&conquer).
     * 2) Probalby a distinction should be made between small matrices (element scale matrices, Jacobians, etc) and
     * large matrices (subdomain scale matrices). It is unreasonable to expect the user to go through a special matrix
     * multiplication algorithm for B^T*E*B, where E is around (6x6), or not be able to inverse a 3x3 Jacobian matrix.
     * In that respect, any sparse formats are not needed for small matrices (unless to save up on RAM if they are 
     * cached). Dense and DenseSymmetric should dominate that category. On the other hand DenseMatrices for small 
     * subdomain matrices are not too far-fetched (Matlab does it).
     * 3) The above kinda destroys any uniformity, but most engineers know when a matrix is small or large. Besides,
     * the transparency inside these 2 categories might actually be increased. Especially small matrices should be 
     * completely transparent since they are going to be passed around alot. Encapsulation is also increased since
     * any knowledge of HPC stuff (e.g. distributed memory formats, 3rd party libraries, GPU routines) does not apply
     * to small matrices. Additionaly, large matrices can (perhaps should not) include solution algorithms inside their
     * class, instead of having separate classes, such as SkylineSolver, etc. 
     * 4) Mutability for large sparse matrices makes sense for changing values within the pattern or even changing the
     * pattern during construction (e.g DOK; still that could be provided by the concrete DOK class). 
     */
    public interface IMatrixView: IReducible //TODO: Needs square matrix, symmetric, spd interfaces (same storage, different methods)
    {
        double this[int row, int col] { get; }
        int Rows { get; }
        int Columns { get; }
        int NumNonZeros { get; } // structural non zeros

        IMatrix DoPointwise(IMatrixView other, Func<double, double, double> binaryOperation);
        IMatrix DoToAllEntries(Func<double, double> unaryOperation);
        IVector ExtractColumn(int col);
        IVector ExtractDiagonal(); // This is only for square matrices
        IVector ExtractRow(int row); //Should there be an IndexOutOfBounds check for the entries in rows, columns? Same for all extract methods.
        IMatrix ExtractSubmatrix(int[] rows, int[] columns); //unsorted indices will also permute the submatrix!
        IVector MultiplyLeft(IVectorView vector);
        IVector MultiplyRight(IVectorView vector);
        IMatrix MultiplyRight(IMatrixView matrix);
        void Print();
        IMatrix Transpose();

        /// <summary>
        /// No extra memory needed, since the transposed view has the same backing data as this instance. If the 
        /// untransposed matrix (this instance) is row major, then iterating the transposed view in column major order 
        /// (e.g. transposedView.MultiplyLeft()) is more efficient.
        /// </summary>
        /// <returns></returns>
        IMatrixView TransposedView(); 
        void Write(string path);
    }
}
