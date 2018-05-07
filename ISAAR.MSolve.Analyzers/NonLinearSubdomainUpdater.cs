using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearSubdomainUpdater : INonLinearSubdomainUpdater
    {
        private readonly Subdomain subdomain;

        public NonLinearSubdomainUpdater(Subdomain subdomain)
        {
            this.subdomain = subdomain;
        }

        public IVector GetRHSFromSolution(IVector solution, IVector dSolution) //TODO leave 
        {
            var forces = new Vector(this.subdomain.TotalDOFs);
            foreach (Element element in this.subdomain.ElementsDictionary.Values)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria
                double[] localSolution = this.subdomain.CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = this.subdomain.CalculateElementNodalDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.MaterialsModified = true;
                double[] f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                this.subdomain.AddLocalVectorToGlobal(element, f, forces.Data);
            }
            return forces;
        }

        public void ResetState()
        {
            foreach (Element element in this.subdomain.ElementsDictionary.Values) element.ElementType.ClearMaterialStresses();
        }

        public void UpdateState()
        {
            foreach (Element element in this.subdomain.ElementsDictionary.Values) element.ElementType.SaveMaterialState();
        }
    }
}
