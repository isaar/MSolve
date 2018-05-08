using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Algorithms.MenkBordas
{
    class MenkBordasVector
    {
        private readonly int numSubdomains;
        private readonly int numEquations; // Not sure if needed

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

        public double DotProduct(MenkBordasVector other)
        {
            double dot = this.Vs * other.Vs;
            for (int i = 0; i < numSubdomains; ++i) dot += this.Ve[i] * other.Ve[i];
            dot += this.Vc * other.Vc; // TODO: avoid this if one or both are 0
            return dot;
        }
    }
}
