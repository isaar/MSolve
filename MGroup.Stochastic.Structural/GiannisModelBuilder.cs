using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using MGroup.Stochastic.Interfaces;

namespace MGroup.Stochastic.Structural
{
    public class GiannisModelBuilder
    {
        public GiannisModelBuilder()
        {
        }

        public Model GetModel(IUncertainParameterRealizer stochasticRealizer, IStochasticDomainMapper domainMapper, int iteration)
        {


            var m = new Model();
            m.SubdomainsDictionary.Add(0, new Subdomain(0));

            m.NodesDictionary.Add(0, new Node(id: 0, x: 0.0, y:  0, z: 0 ));
            m.NodesDictionary.Add(1, new Node(id: 1, x: 0.1, y:  0, z: 0 ));
            m.NodesDictionary.Add(2, new Node(id: 2, x: 0.2, y:  0, z: 0 ));
            m.NodesDictionary.Add(3, new Node(id: 3, x: 0.3, y:  0, z: 0 ));
            m.NodesDictionary.Add(4, new Node(id: 4, x: 0.4, y:  0, z: 0 ));
            m.NodesDictionary.Add(5, new Node(id: 5, x: 0.5, y:  0, z: 0 ));
            m.NodesDictionary.Add(6, new Node(id: 6, x: 0.6, y:  0, z: 0 ));
            m.NodesDictionary.Add(7, new Node(id: 7, x: 0.7, y:  0, z: 0 ));
            m.NodesDictionary.Add(8, new Node(id: 8, x: 0.8, y:  0, z: 0 ));
            m.NodesDictionary.Add(9, new Node(id: 9, x: 0.9, y:  0, z: 0 ));
            m.NodesDictionary.Add(10, new Node(id: 10, x: 1, y:  0, z: 0 ));

            for (int i = 0; i < m.NodesDictionary.Count - 1; i++)
            {
                var e = new Element()
                {
                    ID = i,
                    ElementType = new EulerBeam3D(stochasticRealizer.Realize(iteration, domainMapper, 
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

            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY, Amount = 0 });
            m.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ, Amount = 0 });

            m.Loads.Add(new Load() { Amount = 10, Node = m.NodesDictionary[10], DOF = StructuralDof.TranslationZ });

            //m.ConnectDataStructures();

            return m;


        }
    }
}
