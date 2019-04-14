using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public class NonLinearSubdomainUpdater : INonLinearSubdomainUpdater
    {
        private readonly ISubdomain subdomain;

        public NonLinearSubdomainUpdater(ISubdomain subdomain)
        {
            this.subdomain = subdomain;
        }

        public void ScaleConstraints(double scalingFactor)
        {
            this.subdomain.ScaleConstraints(scalingFactor);
        }

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            return subdomain.GetRhsFromSolution(solution, dSolution);
        }

        public void ResetState()
        {
            this.subdomain.ClearMaterialStresses();
        }

        public void UpdateState()
        {
            this.subdomain.SaveMaterialState();
        }
    }
}

