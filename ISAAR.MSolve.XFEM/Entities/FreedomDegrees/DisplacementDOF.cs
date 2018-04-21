using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: implement a dynamic ordinal numbering that will change depending what dofs are used .
namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    public class DisplacementDOF: IDOF
    {
        public static readonly DisplacementDOF X = new DisplacementDOF("Displacement along X axis");
        public static readonly DisplacementDOF Y = new DisplacementDOF("Displacement along Y axis");
        public static readonly DisplacementDOF Z = new DisplacementDOF("Displacement along Z axis");

        private readonly string name; 

        private DisplacementDOF(string name) // No more displacement dofs can be created
        {
            this.name = name;
        }

        public override string ToString()
        {
            return name;
        }
    }
}
