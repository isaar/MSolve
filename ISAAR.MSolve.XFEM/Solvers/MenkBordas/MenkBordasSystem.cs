using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

//TODO: this only works for cholesky preconditioning of Kss. Extend it to other precoditioners. To abstract it, use a custom
//      preconditioned matrix class (IDisposable) that does multiplications with Ps, Ps^T and KssBar = Ps^T*Kss*Ps
//TODO: ok, this class might be necessarily coupled with XSubdomain2D to keep track of the stored ones. However,
//      MenkBordasPrecMatrix and MenkBordasPreconditioner could just operate on arrays. Even this one could, using the subdomain
//      IDs and Dictionararies. Granted, IDs are easy to overlap, if the subdomains get updated.
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    //TODO: this class should manage all matrices and vectors while they update. 
    class MenkBordasSystem //TODO: IDisposable. Also dispose of modified matrices and factorizations.
    {
        public MenkBordasSystem(DOKSymmetricColMajor Kss, Vector bs)
        {
            this.numDofsStd = Kss.NumRows;
            this.bs = bs;
            this.Ps = MenkBordasPreconditioner.CreateStandardPreconditioner(Kss);
            Kss.Clear(); // No longer needed.

            this.subdomains = new SortedSet<XSubdomain2D>();
            this.modifiedSubdomains = new Dictionary<XSubdomain2D, bool>();
            this.Kee = new Dictionary<XSubdomain2D, DOKSymmetricColMajor>();
            this.Kes = new Dictionary<XSubdomain2D, CSRMatrix>();
            this.Kse = new Dictionary<XSubdomain2D, CSRMatrix>();
            this.B = new Dictionary<XSubdomain2D, SignedBooleanMatrix>();
            this.be = new Dictionary<XSubdomain2D, Vector>();
            this.Pe = new Dictionary<XSubdomain2D, CholeskySuiteSparse>();
        }

        private readonly int numDofsStd;
        private readonly CholeskySuiteSparse Ps; 
        private readonly Vector bs;
        private readonly SortedSet<XSubdomain2D> subdomains;
        private readonly Dictionary<XSubdomain2D, bool> modifiedSubdomains;
        private readonly Dictionary<XSubdomain2D, DOKSymmetricColMajor> Kee; //TODO: only save the factorizations
        private readonly Dictionary<XSubdomain2D, CSRMatrix> Kes;
        private readonly Dictionary<XSubdomain2D, CSRMatrix> Kse;
        private Dictionary<XSubdomain2D, SignedBooleanMatrix> B;
        private readonly Dictionary<XSubdomain2D, Vector> be;
        private readonly Dictionary<XSubdomain2D, CholeskySuiteSparse> Pe;
        // TODO: investigate if some submatrices of L, Q can be cached.

        public void SetSubdomainMatrices(XSubdomain2D subdomain, DOKSymmetricColMajor Kee, CSRMatrix Kes, CSRMatrix Kse,
            Vector be) 
        {
            // Update tracked subdomains
            bool isNew = subdomains.Add(subdomain);
            modifiedSubdomains[subdomain] = true;

            // Dispose unmanaged factorization data of updated matrices
            if (!isNew) this.Pe[subdomain].Dispose();
            this.Pe.Remove(subdomain);

            // Set the new matrices. The old ones will be GCed.
            this.Kee.Add(subdomain, Kee); // This should haven been removed or not added in the first place.
            this.Kes[subdomain] = Kes;
            this.Kse[subdomain] = Kse;
            this.be[subdomain] = be;
        }

        public void SetBooleanMatrices(IDictionary<XSubdomain2D, SignedBooleanMatrix> B) // TODO: allow some to not change
        {
            this.B = new Dictionary<XSubdomain2D, SignedBooleanMatrix>(B);
        }

        public void RemoveSubdomain(XSubdomain2D subdomain) // Not sure of needed. Once enriched, always enriched.
        {
            // Update tracked subdomains
            bool exists = subdomains.Remove(subdomain);
            if (!exists) throw new KeyNotFoundException("No such subdomain is stored");
            modifiedSubdomains.Remove(subdomain);

            // Dispose unmanaged factorization data of updated matrices
            Pe[subdomain].Dispose();

            // Remove the matrices.
            Kee.Remove(subdomain);
            Kes.Remove(subdomain);
            Kse.Remove(subdomain);
            be.Remove(subdomain);
            Pe.Remove(subdomain);
        }

        public void CheckDimensions() //TODO: update this
        {
            // Standard dofs
            Preconditions.CheckSystemSolutionDimensions(numDofsStd, numDofsStd, bs.Length);

            // Number of subdomains
            int numSubdomains = subdomains.Count;
            if ((Kee.Count != numSubdomains) || (Kes.Count != numSubdomains) || (Kse.Count != numSubdomains) 
                || (B.Count != numSubdomains) || (be.Count != numSubdomains))
            {
                throw new ArgumentException($"Mismatch: {numSubdomains} subdomains were declared, but there are"
                    + $" {Kee.Count} enriched-enriched stiffness matrices, {Kes.Count} enriched-standard stiffness matrices,"
                    + $" {Kse.Count} standard-enriched stiffness matrices, {B.Count} signed boolean matrices and"
                    + $" {be.Count} right hand side vectors.");
            }

            // Enriched dofs
            foreach (var sub in modifiedSubdomains.Keys)
            {
                //TODO: All dimensions could be checked.
                Preconditions.CheckSystemSolutionDimensions(Kee[sub], be[sub]);
                Preconditions.CheckSystemSolutionDimensions(Kes[sub], be[sub]);
                Preconditions.CheckSystemSolutionDimensions(Kse[sub], bs);
                Preconditions.CheckMultiplicationDimensions(B[sub].NumColumns, be[sub].Length); //rhs.Length = lhs.Length
            }
        }

        //public (MenkBordasMatrix matrix, Vector rhs) BuildSystem(DOKSymmetricColMajor Kss)
        //{
        //    // Count various dimensions
        //    Dimensions dim = CountDimensions();

        //    // Create the matrix
        //    var copyKss = DOKRowMajor.CreateFromSparseMatrix(Kss).BuildCSRMatrix(true);
        //    var copyKee = new CSRMatrix[dim.NumSubdomains];
        //    int i = 0;
        //    foreach (var sub in modifiedSubdomains.Keys)
        //    {
        //        copyKee[i++] = DOKRowMajor.CreateFromSparseMatrix(Kee[sub]).BuildCSRMatrix(true);
        //    }
        //    var matrix = new MenkBordasMatrix(dim, copyKss, copyKee, Kes.Values.ToArray(), Kse.Values.ToArray(), B.Values.ToArray());

        //    // Assemble the rhs vector
        //    var rhs = Vector.CreateZero(dim.NumDofsAll);
        //    rhs.CopyFromVector(0, bs, 0, dim.NumDofsStd);
        //    for (i = 0; i < dim.NumSubdomains; ++i)
        //    {
        //        rhs.CopyFromVector(dim.SubdomainStarts[i], be[i], 0, be[i].Length);
        //    }

        //    return (matrix, rhs);
        //}

        public (MenkBordasPrecondMatrix matrix, Vector rhs) BuildPreconditionedSystem() // TODO: cache the matrices
        {
            Dimensions dim = CountDimensions();

            // Copy to arrays
            //Vector[] beArray = be.Values.ToArray();
            //CSRMatrix[] KesArray = Kes.Values.ToArray();
            //CSRMatrix[] KseArray = Kse.Values.ToArray();
            //SignedBooleanMatrix[] BArray = B.Values.ToArray();

            // Assemble the rhs vector // TODO: just set the subvectors that change
            var rhs = Vector.CreateZero(dim.NumDofsAll);
            rhs.CopyFromVector(0, bs, 0, dim.NumDofsStd);
            foreach (var sub in subdomains)
            {
                rhs.SetSubvector(be[sub], dim.SubdomainStarts[sub]);
            }

            // Create the preconditioned enriched matrices
            foreach (var subdomain in modifiedSubdomains)
            {
                if (subdomain.Value) // Do not recreate the preconditioners of unmodified Kee
                {
                    DOKSymmetricColMajor kee = Kee[subdomain.Key];
                    // New subdomain: there in no Pe. Modified subdomain: Pe was disposed & removed in the setter.
                    Pe.Add(subdomain.Key, MenkBordasPreconditioner.CreateEnrichedPreconditioner(kee)); 
                    // Dispose of each Kee, once it is no longer needed.
                    Kee.Remove(subdomain.Key);
                    kee.Clear();
                }
            }
            CholeskySuiteSparse[] PeArray = Pe.Values.ToArray();

            // Handle L,Q matrices
            (Matrix L, Matrix Q) = MenkBordasPreconditioner.CreateContinuityEquationsPreconditioners(dim, B, Pe);
            //TODO: I should probably dispose of B here.

            var matrix = new MenkBordasPrecondMatrix(dim, Kes, Kse, Ps, Pe, L, Q);
            return (matrix, rhs);
        }

        public Dimensions CountDimensions()
        {
            var subdomainStarts = new SortedDictionary<XSubdomain2D, int>();
            var subdomainEnds = new SortedDictionary<XSubdomain2D, int>();
            int nextStart = numDofsStd;
            foreach (var subdomainKee in Kee)
            {
                subdomainStarts.Add(subdomainKee.Key, nextStart);
                nextStart = nextStart + subdomainKee.Value.NumRows;
                subdomainEnds.Add(subdomainKee.Key, nextStart);
            }
            int equationsStart = nextStart;
            int numEquations = 0;
            foreach (var subdomainMatrix in B)
            {
                numEquations = subdomainMatrix.Value.NumRows;
                break;
            }
            int numDofsAll = nextStart + numEquations;
            int numDofsEnr = numDofsAll - numDofsStd - numEquations;

            return new Dimensions(numEquations, numDofsStd, numDofsEnr, numDofsAll, subdomains,
                subdomainStarts, subdomainEnds, equationsStart);
        }

        public class Dimensions
        {
            public Dimensions(int numEquations, int numDofsStd, int numDofsEnr, int numDofsAll,
                SortedSet<XSubdomain2D> subdomains,  IReadOnlyDictionary<XSubdomain2D, int> subdomainStarts, 
                IReadOnlyDictionary<XSubdomain2D, int> subdomainEnds, int equationsStart)
            {
                this.NumSubdomains = subdomains.Count;
                this.NumEquations = numEquations;
                this.NumDofsStd = numDofsStd;
                this.NumDofsEnr = numDofsEnr;
                this.NumDofsAll = numDofsAll;
                this.Subdomains = subdomains;
                this.SubdomainStarts = subdomainStarts;
                this.SubdomainEnds = subdomainEnds;
                this.EquationsStart = equationsStart;
            }

            public int NumSubdomains { get; }
            public int NumEquations { get; }
            public int NumDofsStd { get; }
            public int NumDofsEnr { get; }
            public int NumDofsAll { get; }
            public SortedSet<XSubdomain2D> Subdomains { get; }
            public IReadOnlyDictionary<XSubdomain2D, int> SubdomainStarts { get; }
            public IReadOnlyDictionary<XSubdomain2D, int> SubdomainEnds { get; }
            public int EquationsStart { get; }
        }
    }
}
