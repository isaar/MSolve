using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Load
    {
        public ControlPoint ControlPoint { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
