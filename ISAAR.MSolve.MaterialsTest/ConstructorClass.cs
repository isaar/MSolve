using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public class ConstructorClass
    {
        public void MakeElements()
        {
            var materialProperty = new AmbroseMaterialProperty(1, 1, 1, 1);
            var element = new ElementWithMaterial(materialProperty);
            var element2 = new ElementWithoutMaterial(1, 1);
        }
    }
}
