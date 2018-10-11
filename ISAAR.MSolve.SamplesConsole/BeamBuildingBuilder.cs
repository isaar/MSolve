using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class BeamBuildingBuilder
    {
        private static Node GetNodeUnderNode(Model model, Node node)
        {
            double distance = Double.MaxValue;
            Node returnNode = null;
            foreach (Node n in model.NodesDictionary.Values)
            {
                if (n.ID == node.ID) continue;
                if (n.X != node.X || n.Z != node.Z) continue;
                double currentDistance = node.Y - n.Y;
                if (currentDistance < 0) continue;
                if (distance > currentDistance)
                {
                    distance = currentDistance;
                    returnNode = n;
                }
            }
            return returnNode;
        }

        private static Node[] GetAdjacentNodes(Model model, Node node)
        {
            Node nodeLeft = null;
            Node nodeRight = null;
            Node nodeTop = null;
            Node nodeBottom = null;
            double distanceLeft = Double.MaxValue;
            double distanceRight = Double.MaxValue;
            double distanceTop = Double.MaxValue;
            double distanceBottom = Double.MaxValue;
            foreach (Node n in model.NodesDictionary.Values)
            {
                if (n.ID == node.ID) continue;
                if (n.Y != node.Y) continue;
                if (n.X == node.X)
                {
                    if (n.Z - node.Z > 0)
                    {
                        if (distanceTop > n.Z - node.Z)
                        {
                            distanceTop = n.Z - node.Z;
                            nodeTop = n;
                        }
                    }
                    else
                    {
                        if (distanceBottom > node.Z - n.Z)
                        {
                            distanceBottom = node.Z - n.Z;
                            nodeBottom = n;
                        }
                    }
                }
                if (n.Z == node.Z)
                {
                    if (n.X - node.X > 0)
                    {
                        if (distanceRight > n.X - node.X)
                        {
                            distanceRight = n.X - node.X;
                            nodeRight = n;
                        }
                    }
                    else
                    {
                        if (distanceLeft > node.X - n.X)
                        {
                            distanceLeft = node.X - n.X;
                            nodeLeft = n;
                        }
                    }
                }
            }
            return new Node[] { nodeRight, nodeTop, nodeLeft, nodeBottom };
        }

        public static void MakeBeamBuilding(Model model, double startX, double startY, double startZ, double beamWidth, double beamHeight,
            int startNodeID, int startElementID, int subdomainID, int floors, bool isInHexaSoil, bool hasPiles)
        {
            const int nodesPerFloor = 18;

            int nodeID = startNodeID;
            // Construct node structure
            if (!isInHexaSoil)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 5; k++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + k * beamWidth, Y = startY, Z = startZ + j * beamWidth });
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                        model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotZ);
                        nodeID++;
                    }
                }
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 2 * beamWidth, Y = startY, Z = startZ + 3 * beamWidth });
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotZ);
                nodeID++;
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 3 * beamWidth, Y = startY, Z = startZ + 3 * beamWidth });
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotZ);
                nodeID++;
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 4 * beamWidth, Y = startY, Z = startZ + 3 * beamWidth });
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotZ);
                nodeID++;
            }

            for (int i = 1; i <= floors; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 5; k++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + k * beamWidth, Y = startY + i * beamHeight, Z = startZ + j * beamWidth });
                        nodeID++;
                    }
                }
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 2 * beamWidth, Y = startY + i * beamHeight, Z = startZ + 3 * beamWidth });
                nodeID++;
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 3 * beamWidth, Y = startY + i * beamHeight, Z = startZ + 3 * beamWidth });
                nodeID++;
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 4 * beamWidth, Y = startY + i * beamHeight, Z = startZ + 3 * beamWidth });
                nodeID++;
            }
            List<Node> groundNodes = new List<Node>();
            int elementID = startElementID;
            Element e;
            int fibers = 400;
            double b = 0.3;
            double h = 0.1;
            double youngModulus = 2.1e5;
            double poissonRatio = 0.35;

            if (isInHexaSoil)
            {
                List<Node> sub1Nodes = new List<Node>();
                List<Node> sub2Nodes = new List<Node>();
                // Get ground and sub nodes
                for (int i = 0; i < nodesPerFloor; i++)
                {
                    groundNodes.Add(GetNodeUnderNode(model, model.NodesDictionary[startNodeID + i]));
                    sub1Nodes.Add(GetNodeUnderNode(model, groundNodes[i]));
                    sub2Nodes.Add(GetNodeUnderNode(model, sub1Nodes[i]));
                }

                if (hasPiles)
                {
                    // Create sub2 piles
                    for (int i = 0; i < nodesPerFloor; i++)
                    {
                        Node[] sub1AdjacentNodes = GetAdjacentNodes(model, sub1Nodes[i]);
                        Node[] sub2AdjacentNodes = GetAdjacentNodes(model, sub2Nodes[i]);
                        e = new Element()
                        {
                            ID = elementID,
                            ElementType = new EulerBeam3D(youngModulus, poissonRatio, sub2AdjacentNodes, sub1AdjacentNodes)
                            {
                                Density = 7.85,
                                SectionArea = b * h,
                                MomentOfInertiaY = b * b * b * h,
                                MomentOfInertiaZ = b * h * h * h,
                            }
                        };
                        e.NodesDictionary.Add(sub2Nodes[i].ID, sub2Nodes[i]);
                        e.NodesDictionary.Add(sub1Nodes[i].ID, sub1Nodes[i]);
                        foreach (Node node in sub2AdjacentNodes) e.NodesDictionary.Add(node.ID, node);
                        foreach (Node node in sub1AdjacentNodes) e.NodesDictionary.Add(node.ID, node);
                        model.ElementsDictionary.Add(e.ID, e);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                        elementID++;
                    }
                    // Create sub1 piles
                    for (int i = 0; i < nodesPerFloor; i++)
                    {
                        Node[] sub1AdjacentNodes = GetAdjacentNodes(model, sub1Nodes[i]);
                        Node[] groundAdjacentNodes = GetAdjacentNodes(model, groundNodes[i]);
                        e = new Element()
                        {
                            ID = elementID,
                            ElementType = new EulerBeam3D(youngModulus, poissonRatio, sub1AdjacentNodes, groundAdjacentNodes)
                            {
                                Density = 7.85,
                                SectionArea = b * h,
                                MomentOfInertiaY = b * b * b * h,
                                MomentOfInertiaZ = b * h * h * h,
                            }
                        };
                        e.NodesDictionary.Add(sub1Nodes[i].ID, sub1Nodes[i]);
                        e.NodesDictionary.Add(groundNodes[i].ID, groundNodes[i]);
                        foreach (Node node in sub1AdjacentNodes) e.NodesDictionary.Add(node.ID, node);
                        foreach (Node node in groundAdjacentNodes) e.NodesDictionary.Add(node.ID, node);
                        model.ElementsDictionary.Add(e.ID, e);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                        elementID++;
                    }
                }
            }
            else
            {
                for (int i = 0; i < nodesPerFloor; i++)
                    groundNodes.Add(model.NodesDictionary[startNodeID + i]);
                startNodeID += nodesPerFloor;
            }

            // Create ground-to-1st-floor beams
            for (int i = 0; i < nodesPerFloor; i++)
            {
                Node[] groundAdjacentNodes = GetAdjacentNodes(model, groundNodes[i]);
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, isInHexaSoil ? groundAdjacentNodes : null, null)
                    {
                        Density = 7.85,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(groundNodes[i].ID, groundNodes[i]);
                e.NodesDictionary.Add(startNodeID + i, model.NodesDictionary[startNodeID + i]);
                if (isInHexaSoil)
                    foreach (Node node in groundAdjacentNodes)
                        e.NodesDictionary.Add(node.ID, node);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
            }

            for (int i = 0; i < floors; i++)
            {
                double dens = i == floors - 1 ? 10 : 7.85;
                // Horizontal elements
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        e = new Element()
                        {
                            ID = elementID,
                            ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                            {
                                Density = dens,
                                SectionArea = b * h,
                                MomentOfInertiaY = b * b * b * h,
                                MomentOfInertiaZ = b * h * h * h,
                            }
                        };
                        e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + j * 5 + k, model.NodesDictionary[startNodeID + i * nodesPerFloor + j * 5 + k]);
                        e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + j * 5 + k + 1, model.NodesDictionary[startNodeID + i * nodesPerFloor + j * 5 + k + 1]);
                        model.ElementsDictionary.Add(e.ID, e);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                        elementID++;
                    }
                }
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                    {
                        Density = dens,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 3 * 5, model.NodesDictionary[startNodeID + i * nodesPerFloor + 3 * 5]);
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 3 * 5 + 1, model.NodesDictionary[startNodeID + i * nodesPerFloor + 3 * 5 + 1]);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                    {
                        Density = dens,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 3 * 5 + 1, model.NodesDictionary[startNodeID + i * nodesPerFloor + 3 * 5 + 1]);
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 3 * 5 + 2, model.NodesDictionary[startNodeID + i * nodesPerFloor + 3 * 5 + 2]);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
                // Vertical elements
                for (int j = 0; j < 5; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        e = new Element()
                        {
                            ID = elementID,
                            ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                            {
                                Density = dens,
                                SectionArea = b * h,
                                MomentOfInertiaY = b * b * b * h,
                                MomentOfInertiaZ = b * h * h * h,
                            }
                        };
                        e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + k * 5 + j, model.NodesDictionary[startNodeID + i * nodesPerFloor + k * 5 + j]);
                        e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + (k + 1) * 5 + j, model.NodesDictionary[startNodeID + i * nodesPerFloor + (k + 1) * 5 + j]);
                        model.ElementsDictionary.Add(e.ID, e);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                        elementID++;
                    }
                }
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                    {
                        Density = dens,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 2, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 2]);
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 5, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 5]);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                    {
                        Density = dens,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 3, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 3]);
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 6, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 6]);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
                e = new Element()
                {
                    ID = elementID,
                    ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                    {
                        Density = dens,
                        SectionArea = b * h,
                        MomentOfInertiaY = b * b * b * h,
                        MomentOfInertiaZ = b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 4, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 4]);
                e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + 2 * 5 + 7, model.NodesDictionary[startNodeID + i * nodesPerFloor + 2 * 5 + 7]);
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                elementID++;
                // Floor-to-floor beams
                if (i == floors - 1) continue;
                for (int j = 0; j < nodesPerFloor; j++)
                {
                    e = new Element()
                    {
                        ID = elementID,
                        ElementType = new EulerBeam3D(youngModulus, poissonRatio, null, null)
                        {
                            Density = dens,
                            SectionArea = b * h,
                            MomentOfInertiaY = b * b * b * h,
                            MomentOfInertiaZ = b * h * h * h,
                        }
                    };
                    e.NodesDictionary.Add(startNodeID + i * nodesPerFloor + j, model.NodesDictionary[startNodeID + i * nodesPerFloor + j]);
                    e.NodesDictionary.Add(startNodeID + (i + 1) * nodesPerFloor + j, model.NodesDictionary[startNodeID + (i + 1) * nodesPerFloor + j]);
                    model.ElementsDictionary.Add(e.ID, e);
                    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);
                    elementID++;
                }
            }
        }

    }
}
