using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
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
            Vector resultVs = this.Vs.Axpy(other.Vs, scalar);
            var resultVe = new Vector[numSubdomains];
            for (int i = 0; i < numSubdomains; ++i) resultVe[i] = this.Ve[i].Axpy(other.Ve[i], scalar);
            Vector resultVc =  this.Vc.Axpy(other.Vc, scalar); // TODO: avoid this if one or both are 0
            return new MenkBordasVector(numSubdomains, numEquations, resultVs, resultVe, resultVc);
        }

        public void AxpyIntoThis(double scalar, MenkBordasVector other)
        {
            this.Vs.AxpyIntoThis(other.Vs, scalar);
            for (int i = 0; i < numSubdomains; ++i) this.Ve[i].AxpyIntoThis(other.Ve[i], scalar);
            this.Vc.AxpyIntoThis(other.Vc, scalar); // TODO: avoid this if one or both are 0
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
            copy.CopySubvectorFrom(0, Vs, 0, numStdDofs);
            int offset = numStdDofs;
            for (int i = 0; i < numSubdomains; ++i)
            {
                copy.CopySubvectorFrom(offset, Ve[i], 0, Ve[i].Length);
                offset += Ve[i].Length;
            }
            copy.CopySubvectorFrom(offset, Vc, 0, numEquations);

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

        public double Norm2()
        {
            // TODO: should I compute the individual norms, add their squares and then take the root again?
            return Math.Sqrt(this.DotProduct(this));
        }

        public void ScaleIntoThis(double scalar)
        {
            this.Vs.ScaleIntoThis(scalar);
            for (int i = 0; i < numSubdomains; ++i) this.Ve[i].ScaleIntoThis(scalar);
            this.Vc.ScaleIntoThis(scalar); // TODO: avoid this if it is 0
        }

        public void WriteToConsole()
        {
            var writer = new FullVectorWriter(false);
            writer.ArrayFormat = Array1DFormat.PlainVertical;
            writer.NumericFormat = new GeneralNumericFormat();

            Console.WriteLine("\nStandard: ");
            writer.WriteToConsole(Vs);
            for (int i = 0; i < numSubdomains; ++i)
            {
                Console.WriteLine($"\nEnriched subdomain {i}: ");
                writer.WriteToConsole(Ve[i]);
            }
            Console.WriteLine("\nContinuity equations: ");
            writer.WriteToConsole(Vc);
        }
    }
}
