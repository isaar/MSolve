using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Analyzers.Multiscale
{
    public class Square2DRve : IReferenceVolumeElement
    {
        public void ApplyBoundaryConditions(IStructuralModel_v2 model)
        {
            // TODO: this is wrong
            model.Nodes[0].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            //model.Nodes[model.Nodes.Count - 1].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
        }

        public double CalculateRveVolume()
        {
            throw new NotImplementedException();
        }
    }
}
