using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
    public class PressureBoundaryCondition : LoadingCondition
    {
        public double Value { get; private set; }
        
        public PressureBoundaryCondition(double pressureValue)
        {
            this.Value = pressureValue;
        }
    }
}
