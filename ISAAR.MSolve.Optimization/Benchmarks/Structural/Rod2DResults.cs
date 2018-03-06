using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.Optimization.Benchmarks.Structural
{
    public class Rod2DResults
    {
        private readonly Subdomain subdomain;
        private readonly ILinearSystem linearSystem;

        public Rod2DResults(Subdomain subdomain, ILinearSystem linearSystem)
        {
            this.subdomain = subdomain;
            this.linearSystem = linearSystem;
        }

        public double AxialRod2DStress(Element element)
        {
            double[] localDisplacements = subdomain.GetLocalVectorFromGlobal(element, linearSystem.Solution);
            Rod2D rod = (Rod2D)element.ElementType;
            return  rod.CalculateAxialStress(element, localDisplacements, null);
        }
    }
}
