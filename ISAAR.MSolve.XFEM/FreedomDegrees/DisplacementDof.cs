using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: implement a dynamic ordinal numbering that will change depending what dofs are used .
namespace ISAAR.MSolve.XFEM.FreedomDegrees
{
    public class DisplacementDof: IDof
    {
        public static readonly DisplacementDof X = new DisplacementDof("Displacement along X axis");
        public static readonly DisplacementDof Y = new DisplacementDof("Displacement along Y axis");
        public static readonly DisplacementDof Z = new DisplacementDof("Displacement along Z axis");

        private readonly string name; 

        private DisplacementDof(string name) // No more displacement dofs can be created
        {
            this.name = name;
        }

        public override string ToString()
        {
            return name;
        }
    }
}
