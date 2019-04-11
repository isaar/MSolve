using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Optimization.Benchmarks.Structural
{
    public class Rod2DResults
    {
        private readonly Subdomain_v2 subdomain;
        private readonly ILinearSystem_v2 linearSystem;

        public Rod2DResults(Subdomain_v2 subdomain, ILinearSystem_v2 linearSystem)
        {
            this.subdomain = subdomain;
            this.linearSystem = linearSystem;
        }

        public double AxialRod2DStress(Element_v2 element)
        {
            double[] localDisplacements = subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(element, linearSystem.Solution);
            var rod = (Rod2D_v2)element.ElementType;
            return  rod.CalculateAxialStress(element, localDisplacements, null);
        }
    }
}
