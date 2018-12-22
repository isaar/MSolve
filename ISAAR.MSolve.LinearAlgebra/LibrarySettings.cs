using System;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Providers.Managed;
using ISAAR.MSolve.LinearAlgebra.Providers.MKL;

//TODO: These should be thread-safe. Update the documentation as well.
//TODO: A different approach is to employ the Abstract Factory pattern. Clients would use factories to create matrices, instead
//      of constructors or static factory methods. These factories would have different implementations: MKL, Managed, CUDA, etc.
//      Advantage: the client could use simultaneously more than one providers. Disadvantages: 1) A lot of methods would need 
//      references to the factory objects (this could be circumvented with singletons/ enum classes) 2) What happens for matrices
//      or vectors that need e.g. both a BLAS and a SuiteSparse provider?
//TODO: The initialization should be done at the beginning of the program to avoid interference with timing.
//TODO: The name dll in the comments of this class and the providers is too windows oriented. Once I have installed the native 
//      libraries for other operating systems, I need to stop calling them dlls.
namespace ISAAR.MSolve.LinearAlgebra
{
    /// <summary>
    /// Represents a linear algebra library and its corresponding set of providers.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public enum LinearAlgebraProviderChoice
    {
        /// <summary>
        /// Providers that call one or more linear algebra libraries written in managed C# code. Those libraries could be 3rd 
        /// party, custom ones or any combination thereof. They should be used if the libraries required by the other providers 
        /// are not installed on the user's system.
        /// </summary>
        Managed,

        /// <summary>
        /// Providers that call the highly optimized Intel Math Kernel Library (native dlls). Note that Intel MKL relies heavily 
        /// on OpenMP (and sometimes TBB), which might not be desirable if the user also wants to fully utilize the available 
        /// CPU cores in their own multi-threaded code.
        /// </summary>
        MKL
    }

    /// <summary>
    /// Allows the user to set global settings, such as whether to use managed C# linear algebra libraries or optimized native
    /// ones. Methods exposed here are not thread-safe yet.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LibrarySettings
    {
        private static LinearAlgebraProviderChoice providers;

        static LibrarySettings()
        {
            LinearAlgebraProviders = LinearAlgebraProviderChoice.Managed;
        }

        /// <summary>
        /// Linear algebra providers are classes that wrap calls to custom or 3rd party libraries for linear algebra operations.
        /// This property allows the user to choose which library will be used. E.g. For good performance MKL providers (native 
        /// dlls) should be chosen, but if the user hasn't installed Intel MKL, then they can always fall back into the managed 
        /// C# providers. The default behaviour is to use managed providers, which are comparatively inefficient. It is strongly 
        /// recommended to use other providers, if the required libraries are indeed installed. 
        /// </summary>
        /// <remarks>
        /// Note that some linear algebra classes work only with a specific library, which is usually stated in their name or 
        /// at least in their description. When using these classes, the user must make sure that the required libraries are
        /// installed on their system, regardless of the choice of providers here.
        /// </remarks>
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