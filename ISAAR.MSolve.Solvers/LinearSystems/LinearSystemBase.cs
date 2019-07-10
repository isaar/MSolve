using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: add state management
//TODO: should this and its implementations be internal? Analyzers and providers should work with ILinearSystem only. 
namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// Base implementation of <see cref="ILinearSystem"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix">The type of the system's matrix.</typeparam>
    /// <typeparam name="TVector">The type of the system's right hand side and solution vectors.</typeparam>
    public abstract class LinearSystemBase<TMatrix, TVector> : ILinearSystem
        where TMatrix : class, IMatrix //TODO: perhaps this should be IMatrixView
        where TVector : class, IVector
    {
        protected const int initialSize = int.MinValue;

        protected LinearSystemBase(ISubdomain subdomain)
        {
            this.Subdomain = subdomain;
        }

        IMatrixView ILinearSystem.Matrix
        {
            get => Matrix;
            set
            {
                if ((value.NumRows != this.Size) || (value.NumColumns != this.Size))
                {
                    throw new NonMatchingDimensionsException("The provided matrix does not match the dimensions or pattern of"
                        + " this linear system. Make sure that it is initialization was delegated to this linear system"
                        + " after the latest dof ordering.");
                }
                foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();
                Matrix = (TMatrix)value;
            }
        }

        public HashSet<ISystemMatrixObserver> MatrixObservers { get; } = new HashSet<ISystemMatrixObserver>();

        IVector ILinearSystem.RhsVector
        {
            get => RhsConcrete;
            set
            {
                if (value.Length != this.Size)
                {
                    throw new NonMatchingDimensionsException("The provided vector does not match the dimensions or pattern of"
                        + " this linear system. Make sure that it is initialization was delegated to this linear system"
                        + " after the latest dof ordering.");
                }
                RhsConcrete = (TVector)value;
            }
        }

        public int Size { get; private set; } = initialSize;

        public ISubdomain Subdomain { get; }


        internal TMatrix Matrix { get; set; }

        public TVector RhsConcrete { get; set; }

        IVectorView ILinearSystem.Solution { get => SolutionConcrete; }
        public TVector SolutionConcrete { get; set; }

        public virtual void Reset()
        {
            foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();

            // Override the method if memory needs to be disposed in a more complicated way.
            RhsConcrete = null;
            SolutionConcrete = null;
            Matrix = null;

            if (Subdomain is IAsymmetricSubdomain asymmetricSubdomain)
            {
                if (asymmetricSubdomain.FreeDofRowOrdering == null)
                {
                    throw new InvalidOperationException("The freedom degrees of a subdomain must"
                        + " be ordered before defining the size of its corresponding linear system.");
                }
                if (asymmetricSubdomain.FreeDofColOrdering == null)
                {
                    throw new InvalidOperationException("The freedom degrees of a subdomain must"
                        + " be ordered before defining the size of its corresponding linear system.");
                }
                Size = asymmetricSubdomain.FreeDofColOrdering.NumFreeDofs;
            }
            else
            {
                if (Subdomain.FreeDofOrdering == null)
                {
                    throw new InvalidOperationException("The freedom degrees of a subdomain must"
                        + " be ordered before defining the size of its corresponding linear system.");
                }
                Size = Subdomain.FreeDofOrdering.NumFreeDofs;
            }
        }

        public abstract IVector CreateZeroVector();
    }
}
