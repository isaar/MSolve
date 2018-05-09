using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Algorithms.MenkBordas
{
    class MenkBordasVector
    {
        public readonly int numSubdomains;
        public readonly int numEquations; // Not sure if needed

        public readonly Vector Vs;  // standard dofs
        public readonly Vector[] Ve; // enriched dofs of each subdomain
        public readonly Vector Vc; // continuity equations. TODO: special provisions if it is null. Otherwise use a ZeroVector: IVector.

        public MenkBordasVector(int numSubdomains, int numEquations, Vector Vs, Vector[] Ve, Vector Vc)
        {
            this.numSubdomains = numSubdomains;
            this.numEquations = numEquations;
            this.Vs = Vs;
            this.Ve = Ve;
            this.Vc = Vc;
        }

        public static MenkBordasVector CreateZeroWithSameDimensions(MenkBordasVector original)
        {
            var Vs = Vector.CreateZero(original.Vs.Length);
            var Ve = new Vector[original.numSubdomains];
            for (int i = 0; i < original.numSubdomains; ++i) Ve[i] = Vector.CreateZero(original.Ve[i].Length);
            var Vc = Vector.CreateZero(original.Vc.Length);
            return new MenkBordasVector(original.numSubdomains, original.numEquations, Vs, Ve, Vc);
        }

        public void AddIntoThis(MenkBordasVector other)
        {
            this.Vs.AddIntoThis(other.Vs);
            for (int i = 0; i < numSubdomains; ++i) this.Ve[i].AddIntoThis(other.Ve[i]);
            this.Vc.AddIntoThis(other.Vc); // TODO: avoid this if one or both are 0
        }

        public MenkBordasVector Axpy(double scalar, MenkBordasVector other)
        {
            Vector resultVs = this.Vs.Axpy(scalar, other.Vs);
            var resultVe = new Vector[numSubdomains];
            for (int i = 0; i < numSubdomains; ++i) resultVe[i] = this.Ve[i].Axpy(scalar, other.Ve[i]);
            Vector resultVc =  this.Vc.Axpy(scalar, other.Vc); // TODO: avoid this if one or both are 0
            return new MenkBordasVector(numSubdomains, numEquations, resultVs, resultVe, resultVc);
        }

        public void AxpyIntoThis(double scalar, MenkBordasVector other)
        {
            this.Vs.AxpyIntoThis(scalar, other.Vs);
            for (int i = 0; i < numSubdomains; ++i) this.Ve[i].AxpyIntoThis(scalar, other.Ve[i]);
            this.Vc.AxpyIntoThis(scalar, other.Vc); // TODO: avoid this if one or both are 0
        }

        public MenkBordasVector Copy()
        {
            var VeCopy = new Vector[numSubdomains];
            for (int i = 0; i < numSubdomains; ++i) VeCopy[i] = Ve[i].Copy();
            return new MenkBordasVector(numSubdomains, numEquations, Vs.Copy(), VeCopy, Vc.Copy()); // TODO: avoid copying Vc if it is 0
        }

        public Vector CopyToDense()
        {
            // Dimensions
            int numStdDofs = Vs.Length;
            int length = numStdDofs + numEquations;
            for (int sub = 0; sub < numSubdomains; ++sub) length += Ve[sub].Length;
            var copy = Vector.CreateZero(length);

            // Copy the subvectors
            copy.CopyFromVector(0, Vs, 0, numStdDofs);
            int offset = numStdDofs;
            for (int i = 0; i < numSubdomains; ++i)
            {
                copy.CopyFromVector(offset, Ve[i], 0, Ve[i].Length);
                offset += Ve[i].Length;
            }
            copy.CopyFromVector(offset, Vc, 0, numEquations);

            return copy;
        }

        public double DotProduct(MenkBordasVector other)
        {
            double dot = this.Vs * other.Vs;
            for (int i = 0; i < numSubdomains; ++i) dot += this.Ve[i] * other.Ve[i];
            dot += this.Vc * other.Vc; // TODO: avoid this if one or both are 0
            return dot;
        }

        public double[] IndividualDots(MenkBordasVector other)
        {
            var dots = new double[numSubdomains + 2];
            dots[0] = this.Vs * other.Vs;
            for (int i = 0; i < numSubdomains; ++i) dots[1 + i] = this.Ve[i] * other.Ve[i];
            dots[numSubdomains + 1] = this.Vc * other.Vc;
            return dots;
        }

        public double[] IndividualNorms2()
        {
            var norms = new double[numSubdomains + 2];
            norms[0] = Vs.Norm2();
            for (int i = 0; i < numSubdomains; ++i) norms[1 + i] = Ve[i].Norm2();
            norms[numSubdomains + 1] = Vc.Norm2();
            return norms;
        }

        public void WriteToConsole()
        {
            var formatting = new Array1DFormatting("\n", "", "\n");
            FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            Console.Write("\nStandard: ");
            (new FullVectorWriter(Vs, false, formatting)).WriteToConsole();
            for (int i = 0; i < numSubdomains; ++i)
            {
                Console.Write($"\nEnriched subdomain {i}: ");

                (new FullVectorWriter(Ve[i], false, formatting)).WriteToConsole();
            }
            Console.Write("\nContinuity equations: ");
            (new FullVectorWriter(Vc, false, formatting)).WriteToConsole();
        }
    }
}
