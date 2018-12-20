using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Providers;

//TODO: These should be thread-safe.
//TODO: The providers might be better as singletons.
//TODO: A different approach is to employ the Abstract Factory pattern. Clients would use factories to create matrices, instead
//      of constructors or static factory methods. These factories would have different implementations: MKL, Managed, CUDA, etc.
//      Advantage: the client could use simultaneously more than one providers. Disadvantages: 1) A lot of methods would need 
//      references to the factory objects (this could be circumvented with singletons/ enum classes) 2) What happens for matrices
//      or vectors that need e.g. both a BLAS and a SuiteSparse provider?
namespace ISAAR.MSolve.LinearAlgebra
{
    public static class LibrarySettings
    {
        public static void SetBlas(ICBlasProvider blasProvider)
        {
            CBlas = blasProvider;
        }

        public static void SetSparseBlas(ISparseBlasProvider sparseBlas)
        {
            SparseBlas = sparseBlas;
        }

        public static void SetLapack(ILapackProvider lapackProvider)
        {
            LapackLinearEquations = new LapackLinearEquationsFacade(lapackProvider);
            LapackLeastSquares = new LapackLeastSquaresFacadeDouble(lapackProvider);
        }

        internal static ICBlasProvider CBlas { get; private set; } = MklCBlasProvider.UniqueInstance;

        internal static ISparseBlasProvider SparseBlas { get; private set; } = MklSparseBlasProvider.UniqueInstance;

        internal static LapackLinearEquationsFacade LapackLinearEquations { get; private set; } 
            = new LapackLinearEquationsFacade(MklLapackProvider.UniqueInstance);

        internal static LapackLeastSquaresFacadeDouble LapackLeastSquares { get; private set; }
            = new LapackLeastSquaresFacadeDouble(MklLapackProvider.UniqueInstance);
    }
}
