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
        public static ICBlasProvider CBlas { get; set; } = ManagedCBlasProvider.UniqueInstance;
        public static ISparseBlasProvider SparseBlas { get; set; } = ManagedSparseBlasProvider.UniqueInstance;
        public static ILapackProvider Lapack { get; set; } = ManagedLapackProvider.UniqueInstance;
    }
}
