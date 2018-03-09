using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
    public class ControlPointLoad : Load
    {
        public ControlPoint ControlPoint { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
