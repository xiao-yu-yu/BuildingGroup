# -*- coding: utf-8 -*-
import os
from qgis.core import *
import math
import numpy as np
from scipy.spatial import ConvexHull
import scipy.spatial
import triangle as tr
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import time

def load_vector_layer(path, layer_name):
    layer = QgsVectorLayer(path, layer_name, 'ogr')
    if not layer.isValid():
        raise ValueError(f"Invalid layer at path: {path}")
    return layer

class OBBOject:
    def __init__(self):
        self.u0, self.u1, self.c = [0, 0], [0, 0], [0, 0]
        self.e0, self.e1 = 0.0, 0.0
        # *_i represents index.
        self.minX_i, self.maxX_i = 0.0, 0.0
        self.minY_i, self.maxY_i = 0.0, 0.0

    # convert the descriptors to rectangle points.
    def toVertexes(self, isClose=True):
        vertexex = []
        vertexex.append([self.c[0] + self.u0[0] * self.e0 + self.u1[0] * self.e1,
                         self.c[1] + self.u0[1] * self.e0 + self.u1[1] * self.e1])
        vertexex.append([self.c[0] + self.u0[0] * self.e0 - self.u1[0] * self.e1,
                         self.c[1] + self.u0[1] * self.e0 - self.u1[1] * self.e1])
        vertexex.append([self.c[0] - self.u0[0] * self.e0 - self.u1[0] * self.e1,
                         self.c[1] - self.u0[1] * self.e0 - self.u1[1] * self.e1])
        vertexex.append([self.c[0] - self.u0[0] * self.e0 + self.u1[0] * self.e1,
                         self.c[1] - self.u0[1] * self.e0 + self.u1[1] * self.e1])
        if isClose:
            vertexex.append([self.c[0] + self.u0[0] * self.e0 + self.u1[0] * self.e1,
                             self.c[1] + self.u0[1] * self.e0 + self.u1[1] * self.e1])
        return vertexex

    # calculate the point (index) that touched rectangle with the maximum X.
    # should pass the cooridates as an argument.
    def pointTouchRectWithMaxX(self, A):
        max_X, max_P = A[self.minX_i][0], self.minX_i
        if A[self.maxX_i][0] > max_X:
            max_X, max_P = A[self.maxX_i][0], self.maxX_i
        if A[self.minY_i][0] > max_X:
            max_X, max_P = A[self.minY_i][0], self.minY_i
        if A[self.maxY_i][0] > max_X:
            max_X, max_P = A[self.maxY_i][0], self.maxY_i
        return max_P

    # def distanceOfPointFromRect(self, P):
    #     vertexex = self.toVertexes()
    #     min_Dis = pointToLine(P, vertexex[3], vertexex[0])
    #     for i in range(0, 3):
    #         if pointToLine(P, vertexex[i], vertexex[i + 1]) < min_Dis:
    #             min_Dis = pointToLine(P, vertexex[i], vertexex[i + 1])
    #     return min_Dis

    def Orientation(self):
        if self.e0 > self.e1:
            if self.u0[1] > 0:
                return math.acos(self.u0[0])
            else:
                return math.pi - math.acos(self.u0[0])
        else:
            if self.u1[1] > 0:
                return math.acos(self.u1[0])
            else:
                return math.pi - math.acos(self.u1[0])

def generate_adjacency_matrix(features, pointtopoint, edgeIDList, partition_feature):
    features_num = len(features)
    adjacency_matrix = np.zeros((features_num, features_num))


    point_indices = np.array([(edgeIDList[p[0]], edgeIDList[p[1]]) for p in pointtopoint])

    adjacency_matrix[point_indices[:, 0], point_indices[:, 1]] = 1
    adjacency_matrix[point_indices[:, 1], point_indices[:, 0]] = 1


    if partition_feature:
        partition_num = [-i for i in range(1, len(partition_feature) + 1)]
        adjacency_matrix = np.delete(adjacency_matrix, partition_num, axis=0)
        adjacency_matrix = np.delete(adjacency_matrix, partition_num, axis=1)

    return adjacency_matrix

def Dot(u, v):
    return u[0] * v[0] + u[1] * v[1]

def center_point(point1, point2):
    return QgsPointXY((point1.x() + point2.x()) / 2, (point1.y() + point2.y()) / 2)

class Point(object):
    def __init__(self,x,y):
        self.x=x
        self.y=y
    def distance(self,otherpoint):
        dx=self.x-otherpoint.x
        dy=self.y-otherpoint.y
        return math.sqrt(dx**2+dy**2)
    def triangle_area(self, p2, p3):
        x1, y1 = self.x,self.y
        x2, y2 = p2.x,p2.y
        x3, y3 = p3.x,p3.y
        return 0.5 * abs(x1 * y2 + x2 * y3 + x3 * y1 - y1 * x2 - y2 * x3 - y3 * x1)

def dett(point1,point2,point3):
    return (point2[0]-point1[0])*(point3[1]-point1[1])-(point2[1]-point1[1])*(point3[0]-point1[0])

def get_mean_radius(A,CX,CY):
    mean_radius=0
    for i in range(0,len(A)-1):
        mean_radius+=math.sqrt(math.pow(A[i][0]-CX,2)+math.pow(A[i][1]-CY,2))
    N=len(A)-1
    return mean_radius/N

def get_basic_parametries_of_Poly(A):
	CX,CY,area,peri=0,0,0,0
	for i in range(0,len(A)-1):
		CX+=A[i][0]
		CY+=A[i][1]
		peri+=math.sqrt(pow(A[i+1][0]-A[i][0],2)+pow(A[i+1][1]-A[i][1],2))
	CX=CX/(len(A)-1)
	CY=CY/(len(A)-1)
	indication_point=A[0]
	for i in range(1,len(A)-1):
		area+=dett(indication_point,A[i],A[i+1])
	return [[CX,CY],abs(area)*0.5,abs(peri)]

def calCN_Junction(node,centerNode,startEndPointEdge,calCenterNodeJunction):
    if node in startEndPointEdge.keys():
        if node not in calCenterNodeJunction.keys():
            calCenterNodeJunction.update({node: [centerNode]})
        else:
            if centerNode not in calCenterNodeJunction[node]:
                calCenterNodeJunction[node].append(centerNode)

def polygon_area(points):
    n = len(points)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += points[i][0] * points[j][1] - points[j][0] * points[i][1]
    return abs(area) / 2

def mininumAreaRectangle(A):
    size, min_area = len(A), 1.7976931348623157e+128
    i, j = 0, size - 1
    d, u0, u1 = [0, 0], [0, 0], [0, 0]
    OBB = OBBOject()
    while i < size:
        length_edge = math.sqrt(pow(A[i][0] - A[j][0], 2) +
                                pow(A[i][1] - A[j][1], 2))
        if length_edge != 0:
            u0[0], u0[1] = (A[i][0] - A[j][0]) / length_edge, (A[i][1] - A[j][1]) / length_edge
            u1[0], u1[1] = 0 - u0[1], u0[0]
            # u0和u1是垂直的且是单位向量,e0,e1是长和宽的一半
            # print(math.sqrt(u0[0]*u0[0]+u0[1]*u0[1]),u0[0]*u1[0]+u0[1]*u1[1])
            min0, max0, min1, max1, minX_i, maxX_i, minY_i, maxY_i = 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0
            for k in range(0, size):
                d[0], d[1] = A[k][0] - A[j][0], A[k][1] - A[j][1]
                # The projection onto the u0,u1.
                dotU0, dotU1 = Dot(d, u0), Dot(d, u1)
                if dotU0 < min0:
                    min0, minX_i = dotU0, k
                if dotU0 > max0:
                    max0, maxX_i = dotU0, k
                if dotU1 < min1:
                    min1, minY_i = dotU1, k
                if dotU1 > max1:
                    max1, maxY_i = dotU1, k

            area = (max0 - min0) * (max1 - min1)
            if area < min_area:
                min_area = area
                # Update the information.
                OBB.c[0] = A[j][0] + (u0[0] * (max0 + min0) + u1[0] * (max1 + min1)) * 0.5
                OBB.c[1] = A[j][1] + (u0[1] * (max0 + min0) + u1[1] * (max1 + min1)) * 0.5
                OBB.u0, OBB.u1 = [u0[0], u0[1]], [u1[0], u1[1]]
                OBB.e0, OBB.e1 = (max0 - min0) * 0.5, (max1 - min1) * 0.5
                OBB.minX_i, OBB.maxX_i, OBB.minY_i, OBB.maxY_i = minX_i, maxX_i, minY_i, maxY_i
        j, i = i, i + 1
    return OBB, min_area

def extracting_geometric_features(building_feature):
    building_geometric_feature, building_boundary = [], []

    for i in range(len(building_feature)):
        building_feature_geometry = building_feature[i].geometry().asMultiPolygon()[0][0]
        building_feature_coordinates = np.array([[point[0], point[1]] for point in building_feature_geometry])
        building_boundary.append(building_feature_coordinates.tolist())

        # Extract basic parameters
        [[CX, CY], area, peri] = get_basic_parametries_of_Poly(building_feature_coordinates)

        # Calculate the normalized coordinates
        building_feature_uniform_coords = building_feature_coordinates - np.array([CX, CY])

        # Calculate circularity and mean radius
        circularity = 4 * math.pi * area / peri ** 2
        mean_radius = get_mean_radius(building_feature_coordinates, CX, CY)

        # Calculate OBB, SMBR, and orientation
        OBB, SMBR_area = mininumAreaRectangle(building_feature_coordinates)
        orientation = OBB.Orientation()

        # Calculate convex hull properties
        convexHull = ConvexHull(building_feature_uniform_coords)
        ConvexHull_peri = convexHull.area
        ConvexHull_area = convexHull.volume

        # # Calculate long edge properties
        # distance = scipy.spatial.distance.cdist(building_feature_uniform_coords[convexHull.vertices],
        #                                         building_feature_uniform_coords[convexHull.vertices], 'euclidean')
        # long_edge_index = np.unravel_index(np.argmax(convexHull.vertices), distance.shape)
        # long_edge_from, long_edge_to = long_edge_index[0], long_edge_index[1]
        # long_edge = np.linalg.norm(
        #     building_feature_uniform_coords[long_edge_to] - building_feature_uniform_coords[long_edge_from])

        # Feature indicating the ratio of the area of the building to the area of its convex hull.
        concavity_area = area / ConvexHull_area

        # Feature indicating the ratio of the perimeter of the building to the perimeter of its convex hull.
        concavity_peri = peri / ConvexHull_peri

        # Feature indicating the ratio of the area of the building to the area of its minimum bounding rectangle (MBR)
        MBRFullness = area / SMBR_area

        # Feature indicating the elongation of the building, computed as the ratio of the longer side of its minimum bounding rectangle to the shorter side.
        Elongation = OBB.e1 / OBB.e0 if OBB.e0 > OBB.e1 else OBB.e0 / OBB.e1

        edge_orien_weight = np.sum(np.arctan2(np.diff(building_feature_uniform_coords[:, 1]),
                                              np.diff(building_feature_uniform_coords[:, 0])) % (
                                               2 * np.pi) * np.linalg.norm(
            np.diff(building_feature_uniform_coords, axis=0), axis=1))

        edge_length_sum = np.sum(np.linalg.norm(np.diff(building_feature_uniform_coords, axis=0), axis=1))

        # Feature indicating the orientation of walls in the building, computed as the weighted average of edge orientations.
        WALL_orien = edge_orien_weight / edge_length_sum if edge_length_sum != 0 else 0

        feature = [area, peri, circularity, mean_radius, orientation, concavity_area, concavity_peri, MBRFullness,
                   Elongation, WALL_orien]
        building_geometric_feature.append(feature)

    return building_geometric_feature, building_boundary

def filt_tris(triangles, pointList, pointList2, edgeIDList):
    revisedTris = []
    polygonList = []
    setPointList2 = set(pointList2)

    # 筛选建筑拓扑错误
    xy = triangles['vertices'].tolist()
    XY = []
    llist = []
    for i in pointList:
        XY.append([i[0], i[1]])
    for i in xy:
        if i not in XY:
            llist.append(i)
    print(llist)

    for tri in triangles['triangles']:
        points = [pointList[i] for i in tri]
        edges = [edgeIDList[i] for i in tri]

        if not setPointList2.issuperset(points) and len(set(edges)) != 1:
            revisedTris.append(tri)

    triangles['triangles'] = np.array(revisedTris)

    for tri in revisedTris:
        polygonList.append(QgsGeometry.fromPolygonXY([[pointList[int(i)] for i in tri]]))

    return revisedTris, len(revisedTris), polygonList

def satics_tris(revisesTris,segmentList,trisNum,pointList):
    nodeToTris = {i: [] for i in range(len(pointList))}
    tris = np.array(revisesTris)
    revisesTrisSets = [set(tri) for tri in revisesTris]

    for tri_index, tri in enumerate(tris):
        for node in tri:
            nodeToTris[node].append(tri_index)

    nerbors = np.full([trisNum, 3], -1)
    segmentSet = {frozenset(i) for i in segmentList}
    pointtopoint = []

    for i in range(trisNum):
        tri = tris[i]
        for j in range(3):
            otherTwoNodes = tri.tolist()
            otherTwoNodes.remove(tri[j])
            pointtopoint.append(otherTwoNodes)
            otherTwoNodesSet = frozenset(otherTwoNodes)

            connectedTriIndices = set()
            for node in otherTwoNodes:
                connectedTriIndices.update(nodeToTris[node])

            for k in connectedTriIndices:
                if k != i and otherTwoNodesSet < revisesTrisSets[k] and otherTwoNodesSet not in segmentSet:
                    nerbors[i, j] = k

    return nerbors, nodeToTris, pointtopoint

def write_shapeFile(input_layer,shapeFileName,polygonList,geometryType = QgsWkbTypes.Polygon):

    fields = QgsFields()
    # crs = QgsProject.instance().crs()
    crs = input_layer.crs()
    transform_context = QgsProject.instance().transformContext()
    save_options = QgsVectorFileWriter.SaveVectorOptions()
    save_options.driverName = "ESRI Shapefile"
    save_options.fileEncoding = "UTF-8"

    writer = QgsVectorFileWriter.create(
        shapeFileName,
        fields,
        geometryType,
        crs,
        transform_context,
        save_options
    )
    if writer.hasError() != QgsVectorFileWriter.NoError:
        print("Error when creating shapefile: ", writer.errorMessage())
    for polygon in polygonList:
        fet = QgsFeature()
        fet.setGeometry(polygon)
        # fet.setAttributes([1, "text"])
        writer.addFeature(fet)
    del writer

def barycenter_point(point1, point2, point3):
    return QgsPointXY((point1.x() + point2.x() + point3.x()) / 3, (point1.y() + point2.y() + point3.y()) / 3)

def data_organize(features):
    segmentList = []
    pointList = []
    edgeIDList = []
    startEndPointEdge = dict()

    lastPointIndex = 0  # Initialize outside the loop to maintain continuity

    for edge_id, feature in enumerate(features):
        pointXYList = feature.geometry().asMultiPolygon()[0][0]

        # Exclude the last point if it's identical to the first (closed loop)
        pointXYList = pointXYList[:-1]

        for i, point in enumerate(pointXYList):
            currentPointIndex = len(pointList)  # Current index in pointList

            # Add the current point and its edge_id
            pointList.append(point)
            edgeIDList.append(edge_id)

            if i > 0:  # Create segments for all but the first point in each feature
                segment = [lastPointIndex, currentPointIndex]
                segmentList.append(segment)

            lastPointIndex = currentPointIndex  # Update lastPointIndex for next iteration

        # Create a segment connecting the last point of this feature to the first point
        if len(pointXYList) > 0:
            segmentList.append([lastPointIndex, len(pointList) - len(pointXYList)])


    vertices = np.array(pointList)
    segments = np.array(segmentList)
    buildDict = dict({'vertices': vertices, 'segments': segments})
    deletKeys = []
    for key, value in startEndPointEdge.items():
        if len(value) == 1:
            deletKeys.append(key)
    for key in deletKeys:
        startEndPointEdge.pop(key)

    return segmentList, pointList, edgeIDList, buildDict

def creat_skeleton(revisesTris,nerbors,pointList,trisNum,edgeIDList,partition_feature):
    polylineList = []
    calCenterNodeJunction = dict()

    mean_dis = {}
    vis_dis = {}

    point_index = {}
    for index, value in enumerate(edgeIDList):
        point_index.setdefault(value, []).append((pointList[index][0], pointList[index][1]))

    if partition_feature is not None:
        maxId = max(edgeIDList)
        partition_num = [maxId - i for i in range(len(partition_feature))]

    for i in range(trisNum):
        tri = revisesTris[i]
        triNerbor = nerbors[i]
        triNodeCoord = [pointList[j] for j in tri]
        nerborTrisToNode = np.array(np.where(triNerbor > -1))
        nerborsNum = nerborTrisToNode.shape[1]

        if partition_feature!=None:
            edgeIDIndex = [pointList.index(i) for i in triNodeCoord]
            edgeID = [edgeIDList[i] for i in edgeIDIndex]
            result = next((elem for elem in partition_num if elem in edgeID), None)
            if result and nerborsNum==3 and len(set(edgeID))==3:
                ID = [0, 1, 2]
                nodeFaceToSegmentID = edgeID.index(result)
                ID.remove(nodeFaceToSegmentID)
                fistNerborSegmentNodeID = ID[0]
                secondNerborSegmentNodeID = ID[1]
                nodeFaceToSegment = triNodeCoord[nodeFaceToSegmentID]
                fistNerborSegmentNode = triNodeCoord[fistNerborSegmentNodeID]
                secondNerborSegmentNode = triNodeCoord[secondNerborSegmentNodeID]
                centerNode = center_point(fistNerborSegmentNode, secondNerborSegmentNode)
                polylineList.append(QgsGeometry.fromPolylineXY([nodeFaceToSegment, centerNode]))
                continue
            elif result:
                continue
        if nerborsNum == 1:
            ID = [0, 1, 2]
            nodeFaceToSegmentID = nerborTrisToNode[0, 0]
            ID.remove(nodeFaceToSegmentID)
            fistNerborSegmentNodeID = ID[0]
            secondNerborSegmentNodeID = ID[1]
            nodeFaceToSegment = triNodeCoord[nodeFaceToSegmentID]
            fistNerborSegmentNode = triNodeCoord[fistNerborSegmentNodeID]
            secondNerborSegmentNode = triNodeCoord[secondNerborSegmentNodeID]
            centerNode = center_point(fistNerborSegmentNode, secondNerborSegmentNode)
            polylineList.append(QgsGeometry.fromPolylineXY([nodeFaceToSegment, centerNode]))
        elif nerborsNum == 2:
            ID = [0, 1, 2]
            firstNodeFaceToSegmentID = nerborTrisToNode[0, 0]
            secondNodeFaceToSegmentID = nerborTrisToNode[0, 1]
            ID.remove(nerborTrisToNode[0, 0])
            ID.remove(nerborTrisToNode[0, 1])
            thirdNodeID = ID[0]
            firstNodeFaceToSegment = triNodeCoord[firstNodeFaceToSegmentID]
            secondNodeFaceToSegment = triNodeCoord[secondNodeFaceToSegmentID]
            thirdNode = triNodeCoord[thirdNodeID]
            P1 = Point(firstNodeFaceToSegment[0], firstNodeFaceToSegment[1])
            P2 = Point(secondNodeFaceToSegment[0], secondNodeFaceToSegment[1])
            P3 = Point(thirdNode[0], thirdNode[1])
            tris = [edgeIDList[m] for m in revisesTris[i]]
            edgelink = []
            for elem in tris:
                if elem not in edgelink:
                    edgelink.append(elem)
            edgelink2 = [edgelink[1], edgelink[0]]
            single_mean_distance = (Point.distance(P1, P3) + Point.distance(P2, P3)) / 2.0
            # single_tri_area = Point.triangle_area(P1, P2, P3)
            centerNode1 = center_point(firstNodeFaceToSegment, thirdNode)
            centerNode2 = center_point(secondNodeFaceToSegment, thirdNode)
            center_node1 = Point(centerNode1[0], centerNode1[1])
            center_node2 = Point(centerNode2[0], centerNode2[1])
            single_vis_distance = Point.distance(center_node1, center_node2)
            if tuple(edgelink) in mean_dis or tuple(edgelink2) in mean_dis:
                mean_dis[tuple(edgelink)].append(single_mean_distance)
                mean_dis[tuple(edgelink2)].append(single_mean_distance)
            else:
                mean_dis[tuple(edgelink)] = [single_mean_distance]
                mean_dis[tuple(edgelink2)] = [single_mean_distance]
            if tuple(edgelink) in vis_dis or tuple(edgelink2) in vis_dis:
                vis_dis[tuple(edgelink)].append(single_vis_distance)
                vis_dis[tuple(edgelink2)].append(single_vis_distance)
            else:
                vis_dis[tuple(edgelink)] = [single_vis_distance]
                vis_dis[tuple(edgelink2)] = [single_vis_distance]
            # if tuple(edgelink) in convex_dis or tuple(edgelink2) in convex_dis:
            #     convex_dis[tuple(edgelink)].append(single_tri_area)
            #     convex_dis[tuple(edgelink2)].append(single_tri_area)
            # else:
            #     convex_dis[tuple(edgelink)] = [single_tri_area]
            #     convex_dis[tuple(edgelink2)] = [single_tri_area]
            polylineList.append(QgsGeometry.fromPolylineXY([centerNode1, centerNode2]))
            # if thirdNode in startEndPointEdge.keys():
            #     if thirdNode not in calCenterNodeJunction.keys():
            #         calCenterNodeJunction.update({thirdNode: [centerNode1, centerNode2]})
            #     else:
            #         if centerNode1 not in calCenterNodeJunction[thirdNode]:
            #             calCenterNodeJunction[thirdNode].append(centerNode1)
            #         if centerNode2 not in calCenterNodeJunction[thirdNode]:
            #             calCenterNodeJunction[thirdNode].append(centerNode2)
            # if firstNodeFaceToSegment in startEndPointEdge.keys():
            #     if firstNodeFaceToSegment not in calCenterNodeJunction.keys():
            #         calCenterNodeJunction.update({firstNodeFaceToSegment: [centerNode1]})
            #     else:
            #         if centerNode1 not in calCenterNodeJunction[firstNodeFaceToSegment]:
            #             calCenterNodeJunction[firstNodeFaceToSegment].append(centerNode1)
            # if secondNodeFaceToSegment in startEndPointEdge.keys():
            #     if secondNodeFaceToSegment not in calCenterNodeJunction.keys():
            #         calCenterNodeJunction.update({secondNodeFaceToSegment: [centerNode2]})
            #     else:
            #         if centerNode2 not in calCenterNodeJunction[secondNodeFaceToSegment]:
            #             calCenterNodeJunction[secondNodeFaceToSegment].append(centerNode2)
            # calCN_Junction(firstNodeFaceToSegment, centerNode1, startEndPointEdge, calCenterNodeJunction)
            # calCN_Junction(secondNodeFaceToSegment, centerNode2, startEndPointEdge, calCenterNodeJunction)
        elif nerborsNum == 3:
            firstNode = triNodeCoord[0]
            secondNode = triNodeCoord[1]
            thirdNode = triNodeCoord[2]
            barycenter = barycenter_point(firstNode, secondNode, thirdNode)
            centerNode1 = center_point(firstNode, secondNode)
            centerNode2 = center_point(firstNode, thirdNode)
            centerNode3 = center_point(thirdNode, secondNode)
            polylineList.append(QgsGeometry.fromPolylineXY([barycenter, centerNode1]))
            polylineList.append(QgsGeometry.fromPolylineXY([barycenter, centerNode2]))
            polylineList.append(QgsGeometry.fromPolylineXY([barycenter, centerNode3]))
            # calCN_Junction(firstNode, centerNode1, startEndPointEdge, calCenterNodeJunction)
            # calCN_Junction(secondNode, centerNode2, startEndPointEdge, calCenterNodeJunction)
            # calCN_Junction(thirdNode, centerNode3, startEndPointEdge, calCenterNodeJunction)
    for key in mean_dis:
        mean_dis[key] = sum(mean_dis[key]) / len(mean_dis[key])
    for key in vis_dis:
        vis_dis[key] = sum(vis_dis[key])
    # for key in convex_dis:
    #     values3 = convex_dis[key]
    #     avg_value3 = sum(values3)
    #     feature1_index,feature2_index=key[0],key[1]
    #     feature1_coord,feature2_coord=point_index[feature1_index],point_index[feature2_index]
    #     feature1_area,feature2_area=polygon_area(feature1_coord),polygon_area(feature2_coord)
    #     convex_coord=feature1_coord+feature2_coord
    #     convex=ConvexHull(convex_coord)
    #     new_value=(feature1_area+feature2_area+avg_value3)/convex.volume
    #     convex_dis[key] = new_value
    # weight_dis = {}
    # for key in vis_dis.keys() | mean_dis.keys():
    #     weight_dis[key] = vis_dis.get(key, 0) / mean_dis.get(key, 0)

    return calCenterNodeJunction,polylineList,mean_dis,vis_dis

def extract_edges_from_adjacency(adjacency_matrix):
    edges = []
    for i in range(len(adjacency_matrix)):
        for j in range(i + 1, len(adjacency_matrix)):
            if adjacency_matrix[i, j] == 1:
                edges.append([i, j])
                edges.append([j, i])
    return edges

def revise_skeleton(segmentList,startEndPointEdge,calCenterNodeJunction,pointList,nodeToTris,revisesTris,nerbors,polylineList):
    segmentArray = np.array(segmentList).reshape(-1)
    _segmentList = np.array(segmentList)
    for key, value in startEndPointEdge.items():
        if len(value) == 0:
            continue
        else:
            processNum = len(value)
            try:
                keyCenNodeList = calCenterNodeJunction[key]
            except:
                print(key)
                continue
            distantListToKey = [QgsGeometry.fromPolylineXY([key, i]).length() for i in keyCenNodeList]
        if processNum == 1:
            minIndex = distantListToKey.index(min(distantListToKey))
            polylineList.append(QgsGeometry.fromPolylineXY([key, keyCenNodeList[minIndex]]))
        else:
            keyNodeID = pointList.index(key)
            connectSegments = [_segmentList[i // 2, :].tolist() for i in
                               np.where(segmentArray == keyNodeID)[0].tolist()]
            connectSegmentsArrar = copy.deepcopy(connectSegments)
            keyNodeIDTris = nodeToTris[keyNodeID]
            connectSegmentsTris = []
            # 查询以该汇流点为顶点的三角形
            connectSegmentsTris = copy.deepcopy(keyNodeIDTris)
            # for connectSegment in connectSegments:
            #     connectSegment.remove(keyNodeID)
            #     otherIDTris = nodeToTris[connectSegment[0]]
            #     connectSegmentsTris = connectSegmentsTris + list(set(keyNodeIDTris).intersection(set(otherIDTris)))
            # connectSegmentsTris = list(set(connectSegmentsTris))
            tri1 = []
            for i in connectSegmentsTris:
                if len(np.where(nerbors[i] > -1)[0]) == 1:
                    if nerbors[i][revisesTris[i].index(keyNodeID)] > -1:
                        tri1.append(revisesTris[i])
            tri1EdgeCount = [0 for i in range(len(connectSegmentsArrar))]
            tri1Edge = []
            for i in range(len(connectSegmentsArrar)):
                cs = connectSegmentsArrar[i]
                for tri in tri1:
                    if set(cs) < set(tri):
                        tri1EdgeCount[i] += 1
                        tri1Edge.append(cs)

            # <editor-fold desc="老版本的寻找targetEdge">
            # targetEdge = []
            #筛选不是1三角形的边
            # for edge in tri1Edge:
            #     if edge not in connectSegmentsArrar:
            #         targetEdge.append(edge)
            # </editor-fold>

            targetEdgeIndex = [i for i in range(len(tri1EdgeCount)) if tri1EdgeCount[i]==0]
            targetEdge = [connectSegmentsArrar[i] for i in targetEdgeIndex]
            # 存在targetEdge和不存在的情况，targetEdge是未参与构建1三角形的边
            if len(targetEdge) == 0:
                targetEdge = [connectSegmentsArrar[i] for i in range(len(connectSegmentsArrar)-1)]
            i = 0
            while i < processNum:
                try:
                    minIndex = distantListToKey.index(min(distantListToKey))
                except:
                    print("minIndex is wrong:{0}".format(key))
                    break
                centernode = keyCenNodeList[minIndex]
                polylineList.append(QgsGeometry.fromPolylineXY([key, centernode]))
                tag=[]
                for et in targetEdge:
                    tag.append(judge_node_segment_direction(centernode, [pointList[et[0]], pointList[et[1]]]))
                # tag = judge_node_segment_direction(centernode, [pointList[targetEdge[0]], pointList[targetEdge[1]]])
                keyCenNodeList.remove(centernode)
                distantListToKey.remove(distantListToKey[minIndex])
                deletNodes = []
                for node in keyCenNodeList:
                    _tag = []
                    for et in targetEdge:
                        _tag.append(judge_node_segment_direction(node, [pointList[et[0]], pointList[et[1]]]))
                    # _tag = judge_node_segment_direction(node, [pointList[targetEdge[0]], pointList[targetEdge[1]]])
                    if tag == _tag:
                        deletNodes.append(node)
                for node in deletNodes:
                    distantListToKey.remove(distantListToKey[keyCenNodeList.index(node)])
                    keyCenNodeList.remove(node)
                i += 1
    return polylineList

def standardize_geometric_features(building_geometric_feature):
    # Extract geometric features
    area, peri, circularity, mean_radius, orientation, concavity_area, concavity_peri, MBRFullness, Elongation, WALL_orien = [], [], [], [], [], [], [], [], [], []
    for i in range(len(building_geometric_feature)):
        area.append(building_geometric_feature[i][0])
        peri.append(building_geometric_feature[i][1])
        circularity.append(building_geometric_feature[i][2])
        mean_radius.append(building_geometric_feature[i][3])
        orientation.append(building_geometric_feature[i][4])
        concavity_area.append(building_geometric_feature[i][5])
        concavity_peri.append(building_geometric_feature[i][6])
        MBRFullness.append(building_geometric_feature[i][7])
        Elongation.append(building_geometric_feature[i][8])
        WALL_orien.append(building_geometric_feature[i][9])

    # Standardize features using StandardScaler
    scaler = StandardScaler()
    area = scaler.fit_transform(np.array(area).reshape(-1, 1))
    peri = scaler.fit_transform(np.array(peri).reshape(-1, 1))
    circularity = scaler.fit_transform(np.array(circularity).reshape(-1, 1))
    mean_radius = scaler.fit_transform(np.array(mean_radius).reshape(-1, 1))
    orientation = scaler.fit_transform(np.array(orientation).reshape(-1, 1))
    concavity_area = scaler.fit_transform(np.array(concavity_area).reshape(-1, 1))
    concavity_peri = scaler.fit_transform(np.array(concavity_peri).reshape(-1, 1))
    MBRFullness = scaler.fit_transform(np.array(MBRFullness).reshape(-1, 1))
    Elongation = scaler.fit_transform(np.array(Elongation).reshape(-1, 1))
    WALL_orien = np.array(WALL_orien).reshape(-1, 1)

    # Combine features into a single array
    x = np.array(
        [area, peri, circularity, mean_radius, orientation, concavity_area, concavity_peri, MBRFullness, Elongation,
         WALL_orien]).T
    x = x.reshape(x.shape[1], x.shape[2])

    return x

def calculate_edge_feature(edge, building_boundary, vis_dis, mean_dis):
    edge_distance = {}
    for i in range(len(edge)):
        from_point, to_point = edge[i][0], edge[i][1]
        distance, min_edge = scipy.spatial.distance.cdist(building_boundary[from_point],
                                                          building_boundary[to_point],
                                                          'euclidean'), 1.7976931348623157e+128
        for i in range(0, len(distance)):
            for j in range(0, len(distance[i])):
                if distance[i, j] < min_edge:
                    min_edge = distance[i, j]
        key = [from_point, to_point]
        key = tuple(key)
        edge_distance[key] = min_edge
    edge_feature = {}

    for key in edge_distance.keys():
        if key in vis_dis.keys()  and key in mean_dis.keys():
            edge_feature[key] = [edge_distance[key], mean_dis[key], vis_dis[key]]


    return [[*key, *value] for key, value in edge_feature.items()]

def building_triangle_network(building_feature,partition_feature,partition_path,building_boundary):

    features = building_feature + partition_feature

    segmentList, pointList, edgeIDList, buildDict = data_organize(features)
    segmentList_partition, pointList_partition, edgeIDList_partition, buildDict_partition = data_organize(partition_feature)

    print("Triangulate and filter triangles...")
    trainagules = tr.triangulate(buildDict, "pc")

    revisesTris, trisNum, polygonList = filt_tris(trainagules, pointList, pointList_partition, edgeIDList)

    print("Statics on triangles...")
    nerbors, nodeToTris, pointtopoint = satics_tris(revisesTris, segmentList, trisNum, pointList)
    adjacency_matrix = generate_adjacency_matrix(features, pointtopoint, edgeIDList, partition_feature)
    edge = extract_edges_from_adjacency(adjacency_matrix)

    print("Generated skeleton line...")
    calCenterNodeJunction, polylineList, mean_dis, vis_dis = creat_skeleton(revisesTris,nerbors,pointList,trisNum,edgeIDList,partition_feature)
    # polylineList = revise_skeleton(segmentList, startEndPointEdge, calCenterNodeJunction, pointList, nodeToTris,revisesTris, nerbors, polylineList)
    for i in range(len(partition_feature)):
        partition_line = partition_feature[i].geometry().asMultiPolygon()[0][0]
        polylineList.append(QgsGeometry.fromPolylineXY(partition_line))

    edge_feature = calculate_edge_feature(edge, building_boundary, vis_dis, mean_dis)

    return edge_feature, polygonList, polylineList

def main(folder):

    print("Open the folder:{}".format(folder))
    print("Data input...")
    building_path=os.path.join(folder,"stuAll")+".shp"
    partition_path=os.path.join(folder,"partition")+".shp"
    tri_path = os.path.join(folder, "tri") + ".shp"
    skeleton_path = os.path.join(folder, "skeleton") + ".shp"

    building_layer = load_vector_layer(building_path, 'building')
    partition_layer = load_vector_layer(partition_path, 'partition')

    building_feature=list(building_layer.dataProvider().getFeatures())
    partition_feature = list(partition_layer.dataProvider().getFeatures())

    print("Number of buildings:{}".format(len(building_feature)))
    print("Number of partitions:{}".format(len(partition_feature)))

    print("Extracting geometric features...")
    building_geometric_feature,building_boundary=extracting_geometric_features(building_feature)

    x=building_geometric_feature

    # print("Feature normalization...")
    # x = standardize_geometric_features(building_geometric_feature)
    # print(x.shape,building_geometric_feature.shape)


    print("Constructing triangulation network...")
    edge_feature,polygonList,polylineList=building_triangle_network(building_feature,partition_feature,partition_path,building_boundary)

    print("Output of calculation results...")
    edge_txt = open("./feature/edge.txt", "w")
    for i in edge_feature:
        edge_txt.write(str(i) + '\n')
    edge_txt.close()
    np.save("./feature/building.npy",x)
    write_shapeFile(building_layer, shapeFileName=tri_path, polygonList=polygonList, geometryType=QgsWkbTypes.Polygon)
    write_shapeFile(building_layer, shapeFileName=skeleton_path, polygonList=polylineList,
                    geometryType=QgsWkbTypes.LineString)




if __name__ == "__main__":
    # Replace with your building aggregate folder path
    folder_path = './stu'
    main(folder_path)