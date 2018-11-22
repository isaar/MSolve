using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearSubdomainUpdater_v2 : INonLinearSubdomainUpdater_v2
    {
        private readonly Subdomain_v2 subdomain;

        public NonLinearSubdomainUpdater_v2(Subdomain_v2 subdomain)
        {
            this.subdomain = subdomain;
        }

        public void ScaleConstraints(double scalingFactor)
        {
            this.subdomain.ScaleConstraints(scalingFactor);
        }

        public IVector GetRHSFromSolution(IVectorView solution, IVectorView dSolution) //TODO leave 
        {
            return subdomain.GetRHSFromSolution(solution, dSolution);
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

