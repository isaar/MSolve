using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

//TODO: this only works for cholesky preconditioning of Kss. Extend it to other precoditioners. To abstract it, use a custom
//      preconditioned matrix class (IDisposable) that does multiplications with Ps, Ps^T and KssBar = Ps^T*Kss*Ps
//TODO: ok, this class might be necessarily coupled with XSubdomain2D to keep track of the stored ones. However,
//      MenkBordasPrecMatrix and MenkBordasPreconditioner could just operate on arrays. Even this one could, using the subdomain
//      IDs and Dictionararies. Granted, IDs are easy to overlap, if the subdomains get updated.
//TODO: weird bug. If I use Dictionary instead of SortedDictionary for the data and I do not call ClearSubdomains(), the result 
//      changes. Why?
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    //TODO: this class should manage all matrices and vectors while they update. 
    class MenkBordasSystem : IDisposable
    {
        public MenkBordasSystem(IStandardPreconditionerBuilder PsBuilder, IEnrichedPreconditioning enrichedPreconditioning)
        {
            this.PsBuilder = PsBuilder;
            this.enrichedPreconditioning = enrichedPreconditioning;

            this.subdomains = new SortedSet<XSubdomain2D>();
            this.modifiedSubdomains = new SortedDictionary<XSubdomain2D, bool>();
            this.Kee = new SortedDictionary<XSubdomain2D, DokSymmetric>();
            this.Kes = new SortedDictionary<XSubdomain2D, CsrMatrix>();
            this.Kse = new SortedDictionary<XSubdomain2D, CsrMatrix>();
            this.be = new SortedDictionary<XSubdomain2D, Vector>();
            this.Pe = new SortedDictionary<XSubdomain2D, CholeskySuiteSparse>();
        }

        private readonly IStandardPreconditionerBuilder PsBuilder;
        private readonly IEnrichedPreconditioning enrichedPreconditioning;

        private int numDofsStd;
        private IStandardPreconditioner Ps; 
        private Vector bs;
        private readonly SortedSet<XSubdomain2D> subdomains;
        private readonly SortedDictionary<XSubdomain2D, bool> modifiedSubdomains;
        private readonly SortedDictionary<XSubdomain2D, DokSymmetric> Kee; //TODO: only save the factorizations
        private readonly SortedDictionary<XSubdomain2D, CsrMatrix> Kes;
        private readonly SortedDictionary<XSubdomain2D, CsrMatrix> Kse;
        private SortedDictionary<XSubdomain2D, SignedBooleanMatrix> B;
        private readonly SortedDictionary<XSubdomain2D, Vector> be;
        private readonly SortedDictionary<XSubdomain2D, CholeskySuiteSparse> Pe;
        // TODO: investigate if some submatrices of L, Q can be cached.

        //public void ClearSubdomains() //TODO: also dispose factorizations
        //{
        //    subdomains.Clear();
        //    modifiedSubdomains.Clear();
        //    //Kee.Clear();
        //    //Kes.Clear();
        //    //Kse.Clear();
        //    B = null;
        //    //be.Clear();
        //    Pe.Clear();
        //}

        public void Dispose()
        {
            if (Ps != null) Ps.Dispose();
            if (Pe != null)
            {
                foreach (CholeskySuiteSparse factor in Pe.Values)
                {
                    if (factor != null) factor.Dispose();
                }
            }
        }

        public void ProcessStandardDofs(Vector bs)
        {
            this.numDofsStd = bs.Length;
            this.bs = bs;
            this.Ps = PsBuilder.Build();
        }

        public void SetBooleanMatrices(IDictionary<XSubdomain2D, SignedBooleanMatrix> B) // TODO: allow some to not change
        {
            //Console.WriteLine("Calling MenkBordasSystem.SetBooleanMatrices(B)");
            this.B = new SortedDictionary<XSubdomain2D, SignedBooleanMatrix>(B);
        }

        public void SetSubdomainMatrices(XSubdomain2D subdomain, DokSymmetric Kee, CsrMatrix Kes, CsrMatrix Kse,
            Vector be) 
        {
            // Update tracked subdomains
            bool isNew = subdomains.Add(subdomain);
            modifiedSubdomains[subdomain] = true;

            // Dispose unmanaged factorization data of updated matrices
            if (!isNew)
            {
                this.Pe[subdomain].Dispose();
                this.Pe.Remove(subdomain);
            }

            // Set the new matrices. The old ones will be GCed.
            this.Kee.Add(subdomain, Kee); // This should haven been removed or not added in the first place.
            this.Kes[subdomain] = Kes;
            this.Kse[subdomain] = Kse;
            this.be[subdomain] = be;
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

        //public void CheckDimensions() //TODO: update this
        //{
        //    // Standard dofs
        //    Preconditions.CheckSystemSolutionDimensions(numDofsStd, numDofsStd, bs.Length);

        //    // Number of subdomains
        //    int numSubdomains = subdomains.Count;
        //    if ((Kee.Count != numSubdomains) || (Kes.Count != numSubdomains) || (Kse.Count != numSubdomains) 
        //        || (B.Count != numSubdomains) || (be.Count != numSubdomains))
        //    {
        //        throw new ArgumentException($"Mismatch: {numSubdomains} subdomains were declared, but there are"
        //            + $" {Kee.Count} enriched-enriched stiffness matrices, {Kes.Count} enriched-standard stiffness matrices,"
        //            + $" {Kse.Count} standard-enriched stiffness matrices, {B.Count} signed boolean matrices and"
        //            + $" {be.Count} right hand side vectors.");
        //    }

        //    // Enriched dofs
        //    foreach (var sub in modifiedSubdomains.Keys)
        //    {
        //        //TODO: All dimensions could be checked.
        //        Preconditions.CheckSystemSolutionDimensions(Kee[sub], be[sub]);
        //        Preconditions.CheckSystemSolutionDimensions(Kes[sub], be[sub]);
        //        Preconditions.CheckSystemSolutionDimensions(Kse[sub], bs);
        //        Preconditions.CheckMultiplicationDimensions(B[sub].NumColumns, be[sub].Length); //rhs.Length = lhs.Length
        //    }
        //}

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
            #region debug
            //Console.WriteLine("Num total subdomains = " + subdomains.Count);
            //Console.WriteLine("Num Kee = " + Kee.Count);
            //Console.WriteLine("Num Kes = " + Kes.Count);
            //Console.WriteLine("Num Kse = " + Kse.Count);
            //Console.WriteLine("Num be = " + be.Count);
            //Console.WriteLine("Num B = " + ((B != null) ? B.Count : 0));
            //Console.WriteLine("Num rows in B = " + dim.NumEquations);
            //foreach (var sub_be in be)
            //{
            //    Console.WriteLine("Subdomain = " + sub_be.Key.ID);
            //    (new FullVectorWriter(sub_be.Value, false, Array1DFormat.PlainVertical)).WriteToConsole();
            //    Console.WriteLine();
            //}
            #endregion

            // Assemble the rhs vector // TODO: just set the subvectors that change
            var rhs = Vector.CreateZero(dim.NumDofsAll);
            rhs.CopyFromVector(0, bs, 0, dim.NumDofsStd);
            foreach (var sub in subdomains)
            {
                rhs.SetSubvector(be[sub], dim.SubdomainStarts[sub]);
            }

            // Create the preconditioned enriched matrices
            //Console.WriteLine("Modified subdomains: ");
            foreach (var subdomain in subdomains)
            {
                if (modifiedSubdomains[subdomain]) // Do not recreate the preconditioners of unmodified Kee
                {
                    DokSymmetric kee = Kee[subdomain];

                    #region debug
                    //Console.Write($"Subdomain {subdomain.ID}: ");
                    //for (int i = 0; i < kee.NumColumns; ++i)
                    //{
                    //    if (kee[i, i] == 0.0)
                    //    {
                    //        Console.WriteLine($"Singular kee at index {i}");
                    //    }
                    //}
                    #endregion

                    // New subdomain: there in no Pe. Modified subdomain: Pe was disposed & removed in the setter.
                    //TODO: if it is not discarded in the setter, then this doesn't throw a KeyExists excpetion. Why?
                    CholeskySuiteSparse subPe = enrichedPreconditioning.CreateEnrichedPreconditioner(kee);
                    //Console.WriteLine("Num non zeros in Kee after factorization = " + subPe.NumNonZeros);
                    Pe.Add(subdomain, subPe);

                    // Dispose of each Kee, once it is no longer needed.
                    Kee.Remove(subdomain);
                    kee.Clear(); //TODO: perhaps this should be done in the object that creates Pe

                    // Mark this subdomain us unmodified until the caller modifies it again
                    modifiedSubdomains[subdomain] = false;
                }
            }
            //Console.WriteLine();

            // Handle L,Q matrices
            IFactorizationLQ LQ = null;
            if (B != null)
            {
                LQ = enrichedPreconditioning.CreateContinuityEquationsPreconditioner(dim, B, Pe);
                B = null; // Clear it to make sure an exception is thrown if the caller forgets to update B.
                //TODO: also clear each boolean matrix
            }

            var matrix = new MenkBordasPrecondMatrix(dim, Kes, Kse, Ps, Pe, LQ);
            return (matrix, rhs);
        }

        public Dimensions CountDimensions()
        {
            // Enriched dofs
            var subdomainStarts = new SortedDictionary<XSubdomain2D, int>();
            var subdomainEnds = new SortedDictionary<XSubdomain2D, int>();
            int nextStart = numDofsStd;
            foreach (var subdomainRhs in be)
            {
                subdomainStarts.Add(subdomainRhs.Key, nextStart);
                nextStart = nextStart + subdomainRhs.Value.Length;
                subdomainEnds.Add(subdomainRhs.Key, nextStart);
            }

            // Continuity equations
            if (B != null)
            {
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
            else //TODO: Setting illegal values will at least throw exceptions if sth goes wrong. Find a better way to handle it. 
            {
                int equationsStart = int.MinValue;
                int numEquations = int.MinValue;
                int numDofsAll = nextStart;
                int numDofsEnr = numDofsAll - numDofsStd;
                return new Dimensions(numEquations, numDofsStd, numDofsEnr, numDofsAll, subdomains,
                    subdomainStarts, subdomainEnds, equationsStart);
            }
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
