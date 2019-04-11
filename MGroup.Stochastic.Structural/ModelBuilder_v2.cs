using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;

namespace MGroup.Stochastic.Structural
{
    public class ModelBuilder_v2
    {
        public Model_v2 GetModel(RandomVariable randomVariable, IStochasticDomainMapper domainMapper, int iteration)
        {
            var m = new Model_v2();
            m.NodesDictionary.Add(0, new Node_v2() { ID = 0, X = 0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(1, new Node_v2() { ID = 1, X = 1, Y = 0, Z = 0 });
            m.ElementsDictionary.Add(1, new Element_v2()
            {
                ID = 1,
                ElementType = new EulerBeam3D_v2(randomVariable.Realize(iteration, domainMapper, null), 0.3)
            });
            m.Loads.Add(new Load_v2() { Amount = 10, DOF = DOFType.X, Node = m.NodesDictionary[1] });

            return m;
        }
    }
}
