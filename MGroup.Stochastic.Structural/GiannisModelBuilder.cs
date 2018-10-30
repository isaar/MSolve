using System.Collections.Generic;
using ISAAR.MSolve.FEM;
using MGroup.Stochastic.Structural.StochasticRealizers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;

namespace MGroup.Stochastic.Structural
{
    public class GiannisModelBuilder
    {
        public Model GetModel(KarhunenLoeveCoefficientsProvider stochasticRealizer, int iteration)
        {
            var m = new Model();
            m.NodesDictionary.Add(0, new Node() { ID = 0, X = 0.0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(1, new Node() { ID = 1, X = 0.1, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node() { ID = 2, X = 0.2, Y = 0, Z = 0 });
            m.NodesDictionary.Add(3, new Node() { ID = 3, X = 0.3, Y = 0, Z = 0 });
            m.NodesDictionary.Add(4, new Node() { ID = 4, X = 0.4, Y = 0, Z = 0 });
            m.NodesDictionary.Add(5, new Node() { ID = 5, X = 0.5, Y = 0, Z = 0 });
            m.NodesDictionary.Add(6, new Node() { ID = 6, X = 0.6, Y = 0, Z = 0 });
            m.NodesDictionary.Add(7, new Node() { ID = 7, X = 0.7, Y = 0, Z = 0 });
            m.NodesDictionary.Add(8, new Node() { ID = 8, X = 0.8, Y = 0, Z = 0 });
            m.NodesDictionary.Add(9, new Node() { ID = 9, X = 0.9, Y = 0, Z = 0 });
            m.NodesDictionary.Add(10, new Node() { ID = 10, X = 1, Y = 0, Z = 0 });
            for (int i = 0; i < m.NodesDictionary.Count - 2; i++)
            {
                m.ElementsDictionary.Add(i, new Element()
                {
                    ID = i,
                    ElementType = new EulerBeam3D(stochasticRealizer.Realize(iteration, i), 0.3)
                });
            }

            m.Loads.Add(new Load() { Amount = 10, DOF = DOFType.Z, Node = m.NodesDictionary[10] });

            return m;


        }
    }
}
