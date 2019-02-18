using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public class NonLinearSubdomainUpdater_v2 : INonLinearSubdomainUpdater_v2
    {
        private readonly ISubdomain_v2 subdomain;

        public NonLinearSubdomainUpdater_v2(ISubdomain_v2 subdomain)
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

