using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public interface IMaterialProperty
    {
        IMaterialState BuildMaterialState(double[] coordinates);
    }
}
