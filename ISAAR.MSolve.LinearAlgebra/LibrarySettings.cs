using System;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Providers.Managed;
using ISAAR.MSolve.LinearAlgebra.Providers.MKL;

//TODO: These should be thread-safe.
//TODO: A different approach is to employ the Abstract Factory pattern. Clients would use factories to create matrices, instead
//      of constructors or static factory methods. These factories would have different implementations: MKL, Managed, CUDA, etc.
//      Advantage: the client could use simultaneously more than one providers. Disadvantages: 1) A lot of methods would need 
//      references to the factory objects (this could be circumvented with singletons/ enum classes) 2) What happens for matrices
//      or vectors that need e.g. both a BLAS and a SuiteSparse provider?
//TODO: The initialization should be done at the beginning of the program to avoid interference with timing.
namespace ISAAR.MSolve.LinearAlgebra
{
    public enum LinearAlgebraProviderChoice
    {
        Managed, MKL
    }

    public static class LibrarySettings
    {
        private static LinearAlgebraProviderChoice providers;

        static LibrarySettings()
        {
            LinearAlgebraProviders = LinearAlgebraProviderChoice.Managed;
        }

        public static LinearAlgebraProviderChoice LinearAlgebraProviders
        {
            get => providers;
            set
            {
                providers = value;
                if (providers == LinearAlgebraProviderChoice.Managed)
                {
                    Blas = ManagedBlasProvider.UniqueInstance;
                    BlasExtensions = ManagedBlasExtensionsProvider.UniqueInstance;
                    SparseBlas = ManagedSparseBlasProvider.UniqueInstance;
                    LapackLinearEquations = new LapackLinearEquationsFacade(ManagedLapackProvider.UniqueInstance);
                    LapackLeastSquares = new LapackLeastSquaresFacadeDouble(ManagedLapackProvider.UniqueInstance);
                }
                else if (providers == LinearAlgebraProviderChoice.MKL)
                {
                    Blas = MklBlasProvider.UniqueInstance;
                    BlasExtensions = MklBlasExtensionsProvider.UniqueInstance;
                    SparseBlas = MklSparseBlasProvider.UniqueInstance;
                    LapackLinearEquations = new LapackLinearEquationsFacade(MklLapackProvider.UniqueInstance);
                    LapackLeastSquares = new LapackLeastSquaresFacadeDouble(MklLapackProvider.UniqueInstance);
                }
                else // This is useful to catch errors when adding new providers.
                {
                    throw new NotImplementedException("For now only managed and MKL providers are available.");
                }
            }
        }

        internal static IBlasProvider Blas { get; private set; }

        internal static IBlasExtensionsProvider BlasExtensions { get; private set; } 

        internal static ISparseBlasProvider SparseBlas { get; private set; }

        internal static LapackLinearEquationsFacade LapackLinearEquations { get; private set; } 

        internal static LapackLeastSquaresFacadeDouble LapackLeastSquares { get; private set; }
    }
}
