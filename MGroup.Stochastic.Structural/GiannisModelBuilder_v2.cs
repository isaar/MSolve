using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using MGroup.Stochastic.Interfaces;

namespace MGroup.Stochastic.Structural
{
    public class GiannisModelBuilder_v2
    {
        public GiannisModelBuilder_v2()
        {
        }

        public Model_v2 GetModel(IUncertainParameterRealizer stochasticRealizer, IStochasticDomainMapper domainMapper, int iteration)
        {


            var m = new Model_v2();
            m.SubdomainsDictionary.Add(0, new Subdomain_v2(0));

            m.NodesDictionary.Add(0, new Node_v2() { ID = 0, X = 0.0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(1, new Node_v2() { ID = 1, X = 0.1, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node_v2() { ID = 2, X = 0.2, Y = 0, Z = 0 });
            m.NodesDictionary.Add(3, new Node_v2() { ID = 3, X = 0.3, Y = 0, Z = 0 });
            m.NodesDictionary.Add(4, new Node_v2() { ID = 4, X = 0.4, Y = 0, Z = 0 });
            m.NodesDictionary.Add(5, new Node_v2() { ID = 5, X = 0.5, Y = 0, Z = 0 });
            m.NodesDictionary.Add(6, new Node_v2() { ID = 6, X = 0.6, Y = 0, Z = 0 });
            m.NodesDictionary.Add(7, new Node_v2() { ID = 7, X = 0.7, Y = 0, Z = 0 });
            m.NodesDictionary.Add(8, new Node_v2() { ID = 8, X = 0.8, Y = 0, Z = 0 });
            m.NodesDictionary.Add(9, new Node_v2() { ID = 9, X = 0.9, Y = 0, Z = 0 });
            m.NodesDictionary.Add(10, new Node_v2() { ID = 10, X = 1, Y = 0, Z = 0 });

            for (int i = 0; i < m.NodesDictionary.Count - 1; i++)
            {
                var e = new Element_v2()
                {
                    ID = i,
                    ElementType = new EulerBeam3D_v2(stochasticRealizer.Realize(iteration, domainMapper, 
                        new []
                        {
                            (m.NodesDictionary[i + 1].X + m.NodesDictionary[i].X)/2,
                            (m.NodesDictionary[i + 1].Y + m.NodesDictionary[i].Y)/2,
                            (m.NodesDictionary[i + 1].Z + m.NodesDictionary[i].Z)/2,
                        }), 0.3)
                    {
                        Density = 7.85,
                        SectionArea = 1,
                        MomentOfInertiaY = 1,
                        MomentOfInertiaZ = 1,
                        MomentOfInertiaPolar = 1,
                    }
                };
                e.AddNodes(new[] { m.NodesDictionary[i], m.NodesDictionary[i + 1] });

                m.ElementsDictionary.Add(i, e);
                m.SubdomainsDictionary[0].Elements.Add(e);
            }

            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.X, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.Y, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.Z, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.RotX, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.RotY, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.RotZ, Amount = 0 });

            m.Loads.Add(new Load_v2() { Amount = 10, Node = m.NodesDictionary[10], DOF = DOFType.Z });

            //m.ConnectDataStructures();

            return m;


        }
    }
}
