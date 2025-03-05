#include <iostream>
#include "parser.h"
#include <cmath>
#include "ppm.h"
using namespace std;

/*typedefs*/
typedef struct {
	int first,second,third;
} Vector3i;

typedef struct {
	float first,second,third;
} Vector3f;
typedef struct {
	Vector3f col0,col1,col2;
} Matrix3x3;

typedef struct {
	int first,second;
} Vector2i;

Vector3f fromParserToUsVector3fConverter(parser::Vec3f input){
	Vector3f result;
	result.first = input.x;
	result.second = input.y;
	result.third = input.z;
	return result;
}


class Vector3fFunctions{
	public:
	static Vector3f add (Vector3f parameter1, Vector3f parameter2){
		Vector3f result ;
		result.first = parameter1.first + parameter2.first;
		result.second = parameter1.second + parameter2.second;
		result.third = parameter1.third + parameter2.third;
		return result;
	}
	static Vector3f substract(Vector3f parameter1, Vector3f parameter2){
		Vector3f result ;
		result.first = parameter1.first - parameter2.first;
		result.second = parameter1.second - parameter2.second;
		result.third = parameter1.third - parameter2.third;
		return result;
	}
	static Vector3f negate(Vector3f parameter){
		Vector3f result;
		result.first = -parameter.first;
		result.second = -parameter.second;
		return result;
	}
	static Vector3f multiply(Vector3f parameter, float scalar){
		Vector3f result;
		result.first = scalar*parameter.first;
		result.second = scalar*parameter.second;
		result.third = scalar*parameter.third;
		return result;
	}
	static Vector3f divide(Vector3f parameter, float scalar){
		return multiply(parameter,1/scalar);
	}
	static float magnitude(Vector3f parameter){
		float squareMagnitude = parameter.first*parameter.first+parameter.second*parameter.second+parameter.third*parameter.third;
		return sqrt(squareMagnitude);
	}
	static Vector3f normalize(Vector3f parameter){
		return  divide(parameter,magnitude(parameter));
	}
	static Vector3f crossProduct(Vector3f parameter1 , Vector3f parameter2){
			Vector3f result;
			result.first = parameter1.second*parameter2.third-parameter1.third*parameter2.second;
			result.second = parameter1.third*parameter2.first-parameter1.first*parameter2.third;
			result.third = parameter1.first*parameter2.second-parameter1.second*parameter2.first;

			return result;
	}
	static float dotProduct(Vector3f parameter1, Vector3f parameter2){
		return parameter1.first*parameter2.first+parameter1.second*parameter2.second+parameter1.third*parameter2.third;
	}
	static Vector3f elementwiseMultiplication(Vector3f first, Vector3f second){
		return {first.first*second.first, first.second*second.second, first.third*second.third};
	}
	static Vector3i toIntVectorByLimit(Vector3f inputVector, int limit){
		Vector3i result;
		result.first = round(inputVector.first);
		result.second = round(inputVector.second);
		result.third = round(inputVector.third);

		result.first = min(result.first,limit);
		result.second = min(result.second,limit);
		result.third = min(result.third,limit);
		return result;
	}
	static Vector3f toFloatVector(Vector3i inputVector){
		Vector3f result;
		result.first = (float)inputVector.first;
		result.second = (float)inputVector.second;
		result.third = (float)inputVector.third;
		return result;
	}
	static float computeDeterminant3x3(Matrix3x3 matrix){
		return 
		matrix.col0.first*(matrix.col1.second*matrix.col2.third-matrix.col1.third*matrix.col2.second)
		-
		matrix.col1.first*(matrix.col0.second*matrix.col2.third-matrix.col2.second*matrix.col0.third)
		+
		matrix.col2.first*(matrix.col0.second*matrix.col1.third-matrix.col1.second*matrix.col0.third)
		;
	}
} evaluator;


parser::Scene scene;
struct Ray{
	Vector3f source;
	Vector3f direction;
};

struct IntersectionReport
{
	Vector3f intersectionPoint;
	Vector3f intersectionSurfaceNormal;
	Ray intersectedRay; 
	float t_parameter;
	bool occurance;
};

struct PointLightSource{
	Vector3f position;
	Vector3f color;
};
PointLightSource *pointLightSources;
int pointLightSourceCount;

class Object3D;
Object3D **objects;
int object_count = 0;

Vector3f computeRayColor(Ray ray, int limit);

class Object3D{
	public:
	virtual IntersectionReport getIntersectionReport(Ray ray){throw;}
	int index;
	
	Vector3f ambientReflectance;
	Vector3f diffuseReflectance;
	Vector3f specularReflectance;
	Vector3f mirrorReflectance;
	int phongExp;
	bool isMirrorLike;

	Vector3f ambientContribition;

	Object3D(Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike){
			
			this->ambientReflectance = ambientReflectance;
			this->diffuseReflectance = diffuseReflectance;
			this->specularReflectance = specularReflectance;
			this->mirrorReflectance = mirrorReflectance;
			this->phongExp = phongExp;
			this->isMirrorLike = isMirrorLike;

			Vector3f ambientLight ;
			ambientLight.first  = scene.ambient_light.x;
			ambientLight.second = scene.ambient_light.y;
			ambientLight.third  = scene.ambient_light.z;
			
			ambientContribition = evaluator.elementwiseMultiplication(ambientLight,ambientReflectance);
	}
	
	virtual Vector3f getVisual(IntersectionReport intersectionReport, int limit){	
		Vector3f shadowRayStartingPoint = evaluator.add(intersectionReport.intersectionPoint,
												evaluator.multiply(intersectionReport.intersectionSurfaceNormal,scene.shadow_ray_epsilon));
		
		//ambient component
		Vector3f result = ambientContribition;
		for(int i = 0 ; i < pointLightSourceCount ; i++){
			PointLightSource currentSource = pointLightSources[i];
			//shadowing
			bool prematureIntersection = false;
			
			for(int j = 0 ; j < object_count; j++){
				if(j == this->index){
					continue;
				}
				Ray computedRay;
				computedRay.source = shadowRayStartingPoint;
				computedRay.direction = evaluator.substract(currentSource.position, intersectionReport.intersectionPoint);
				IntersectionReport computedIntersection = objects[j]->getIntersectionReport(computedRay);
				if(computedIntersection.occurance){
					if(computedIntersection.t_parameter < 1){
						prematureIntersection = true;
						break;
					}
				}
			}
			if(prematureIntersection){
				continue;
			}	
			
			Vector3f toCurrentSourceConnector    = evaluator.substract(currentSource.position,intersectionReport.intersectionPoint);
			float toCurrentSourceConnectorMagnitude  = evaluator.magnitude(toCurrentSourceConnector) ;
			Vector3f currentSourceEffectiveColor = evaluator.divide(currentSource.color,
													toCurrentSourceConnectorMagnitude*toCurrentSourceConnectorMagnitude);
				
			Vector3f illuminationVector = evaluator.normalize(toCurrentSourceConnector);
			//diffuse component
			float diffusionCoefficient = evaluator.dotProduct(illuminationVector,intersectionReport.intersectionSurfaceNormal);
			if(diffusionCoefficient > 0){
				Vector3f contribution = evaluator.elementwiseMultiplication(evaluator.multiply(currentSourceEffectiveColor,diffusionCoefficient),diffuseReflectance);
				result = evaluator.add(result,contribution);
			}
			//specular component
			Vector3f to_observer = evaluator.normalize(evaluator.substract(intersectionReport.intersectedRay.source,intersectionReport.intersectionPoint));
			Vector3f halfVector = evaluator.normalize(evaluator.add(to_observer,illuminationVector));
			float temp = evaluator.dotProduct(halfVector, intersectionReport.intersectionSurfaceNormal);
			if(temp > 0){
				float specularContributionCoefficient = pow(temp, phongExp);
				Vector3f specularContribution = evaluator.elementwiseMultiplication(
					specularReflectance,
					evaluator.multiply(currentSourceEffectiveColor, specularContributionCoefficient)
				);
				result = evaluator.add(result,specularContribution);
			}

		}
		//computation of reflectance
		if((limit > 0) && (this->mirrorReflectance.first !=0 || this->mirrorReflectance.second !=0 || this->mirrorReflectance.third !=0) && this->isMirrorLike ){
			Ray computationRequest;
			computationRequest.source = shadowRayStartingPoint;
			computationRequest.direction = evaluator.add(intersectionReport.intersectedRay.direction, 
				evaluator.multiply(intersectionReport.intersectionSurfaceNormal,
					(-2)*evaluator.dotProduct(intersectionReport.intersectedRay.direction,intersectionReport.intersectionSurfaceNormal)
				)
			);
			Vector3f computedRayColor = computeRayColor(computationRequest, limit-1);
			Vector3f computedRayColorContribution = evaluator.elementwiseMultiplication(computedRayColor,this->mirrorReflectance);
			result = evaluator.add(result, computedRayColorContribution);
		}
		return result;
		
	}
};

class Sphere :public Object3D{
	public:
	float radius;
	Vector3f center;

	Sphere(Vector3f parameter_center, 
			float parameter_radius, 
			Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike): 
				Object3D(ambientReflectance, 
				diffuseReflectance,
				specularReflectance, 
				mirrorReflectance ,
				phongExp,
				isMirrorLike){
						radius = parameter_radius;
						center = parameter_center;
	}
	virtual IntersectionReport getIntersectionReport(Ray ray) override{
		
		IntersectionReport result;
		
		float a = evaluator.dotProduct(ray.direction,ray.direction);
		float b = 2*evaluator.dotProduct(ray.direction,evaluator.substract(ray.source,center));
		float c = evaluator.dotProduct(ray.source,ray.source)+evaluator.dotProduct(center,center)-2*evaluator.dotProduct(ray.source,center)-radius*radius;
		float determinant = b*b-4*a*c;
		if(determinant <= 0){
			result.occurance=false;
			return result;
		}
		float sqrt_determinant = sqrt(determinant);
		float t_root1 = (-b+sqrt_determinant)/(2*a);
		float t_root2 = (-b-sqrt_determinant)/(2*a);
		
		if(t_root1 > 0 && t_root2 > 0){
			result.t_parameter = min(t_root1,t_root2);
		}else if(t_root1 > 0 && t_root2 < 0){
			result.t_parameter = t_root1;
		}else if(t_root1 < 0 && t_root2 > 0){
			result.t_parameter = t_root2;
		}else{
			result.occurance = false;
			return result;
		}
		result.occurance = true;
		result.intersectionPoint = evaluator.add(ray.source,evaluator.multiply(ray.direction,result.t_parameter));
		result.intersectionSurfaceNormal = evaluator.normalize(evaluator.substract(result.intersectionPoint,center));
		result.intersectedRay = ray;
		return result;
	}
};

float triangleIntersectionEpsilon = 0.001;

class Triangle: public Object3D{
	private:
		Vector3f v3_to_v1;
		Vector3f v3_to_v2;

		Vector3f v1_to_v2;
		Vector3f v1_to_v3;
		Vector3f normal;
	public:
	Vector3f v1,v2,v3;
	Triangle(Vector3f v1_parameter, Vector3f v2_parameter, Vector3f v3_parameter, 
	
			Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike
	): Object3D(
					ambientReflectance, 
					diffuseReflectance,
					specularReflectance, 
					mirrorReflectance ,
					phongExp,
					isMirrorLike
				){
						v1 = v1_parameter;
						v2 = v2_parameter;
						v3 = v3_parameter;
								
						v3_to_v1 = evaluator.substract(v1,v3);
						v3_to_v2 = evaluator.substract(v2,v3);;

						v1_to_v2 = evaluator.substract(v2,v1);
						v1_to_v3 = evaluator.substract(v3,v1);

						normal = evaluator.normalize(
												evaluator.crossProduct(
													evaluator.substract(v2,v1),
													evaluator.substract(v3,v1)
												)
											);	
				}
	virtual IntersectionReport getIntersectionReport(Ray ray) override{
		IntersectionReport result;
		result.occurance = false;
		
		float c1,c2,t;
		Matrix3x3 matrix = {
			v3_to_v1,
			v3_to_v2,
			evaluator.multiply(ray.direction,-1)
		};
		float determinant = evaluator.computeDeterminant3x3(matrix);
		
		if(determinant == 0){
			return result;
		}else{
			

			Vector3f resultVector = evaluator.substract(ray.source,v3);
			
			Matrix3x3 activeMatrixC1 = matrix;
			activeMatrixC1.col0 = resultVector;
			c1 = evaluator.computeDeterminant3x3(activeMatrixC1)/determinant;
			if(c1 < 0-triangleIntersectionEpsilon || c1 > 1+triangleIntersectionEpsilon ){
				return result;
			}
			Matrix3x3 activeMatrixC2 = matrix;
			activeMatrixC2.col1 = resultVector;
			c2 = evaluator.computeDeterminant3x3(activeMatrixC2)/determinant;
			
			if(c2 < 0-triangleIntersectionEpsilon || c2 > 1+triangleIntersectionEpsilon){
				return result;
			}
			if(1-c1-c2 < 0-triangleIntersectionEpsilon || 1-c1-c2 > 1+triangleIntersectionEpsilon){
				return result;
			}
			
			
			Matrix3x3 activeMatrixT = matrix;
			activeMatrixT.col2 = resultVector;
			t = evaluator.computeDeterminant3x3(activeMatrixT)/determinant;
			if(t < 0){
				return result;
			}
		}
		
		//occurace is guaranteed
		result.occurance = true;
		result.intersectedRay = ray;
		result.intersectionSurfaceNormal = this->normal;

		
		result.t_parameter = t; 
		result.intersectionPoint = evaluator.add(ray.source,
			evaluator.multiply(ray.direction,t)
		);
		return result;
		
	}
};

class AbstractTriangle{
	private:
		Vector3f v3_to_v1;
		Vector3f v3_to_v2;

		Vector3f v1_to_v2;
		Vector3f v1_to_v3;
		Vector3f normal;
		Vector3f v1,v2,v3;
	public:
		AbstractTriangle(){}
		void setup(Vector3f v1_parameter, Vector3f v2_parameter, Vector3f v3_parameter){
					v1 = v1_parameter;
						v2 = v2_parameter;
						v3 = v3_parameter;
								
						v3_to_v1 = evaluator.substract(v1,v3);
						v3_to_v2 = evaluator.substract(v2,v3);;

						v1_to_v2 = evaluator.substract(v2,v1);
						v1_to_v3 = evaluator.substract(v3,v1);

						normal = evaluator.normalize(
												evaluator.crossProduct(
													evaluator.substract(v2,v1),
													evaluator.substract(v3,v1)
												)
											);	
		}
		
		IntersectionReport getIntersectionReport(Ray ray){
		IntersectionReport result;
		result.occurance = false;
		
		float c1,c2,t;
		Matrix3x3 matrix = {
			v3_to_v1,
			v3_to_v2,
			evaluator.multiply(ray.direction,-1)
		};
		float determinant = evaluator.computeDeterminant3x3(matrix);
		
		if(determinant == 0){
			return result;
		}else{
			

			Vector3f resultVector = evaluator.substract(ray.source,v3);
			
			Matrix3x3 activeMatrixC1 = matrix;
			activeMatrixC1.col0 = resultVector;
			c1 = evaluator.computeDeterminant3x3(activeMatrixC1)/determinant;
			if(c1 < 0-triangleIntersectionEpsilon || c1 > 1+triangleIntersectionEpsilon ){
				return result;
			}
			Matrix3x3 activeMatrixC2 = matrix;
			activeMatrixC2.col1 = resultVector;
			c2 = evaluator.computeDeterminant3x3(activeMatrixC2)/determinant;
			
			if(c2 < 0-triangleIntersectionEpsilon || c2 > 1+triangleIntersectionEpsilon){
				return result;
			}
			if(1-c1-c2 < 0-triangleIntersectionEpsilon || 1-c1-c2 > 1+triangleIntersectionEpsilon){
				return result;
			}
			
			
			Matrix3x3 activeMatrixT = matrix;
			activeMatrixT.col2 = resultVector;
			t = evaluator.computeDeterminant3x3(activeMatrixT)/determinant;
			if(t < 0){
				return result;
			}
		}
		
		//occurace is guaranteed
		result.occurance = true;
		result.intersectedRay = ray;
		result.intersectionSurfaceNormal = this->normal;

		
		result.t_parameter = t; 
		result.intersectionPoint = evaluator.add(ray.source,
			evaluator.multiply(ray.direction,t)
		);
		return result;
		
	}
		

};

Vector3f getVertex(int parsed_type_id){
	Vector3f result ;
	result = fromParserToUsVector3fConverter(scene.vertex_data[parsed_type_id-1]);
	return result;
}
class Mesh3D: public Object3D{
	private:
	Mesh3D **children;
	bool hasChildren;

	Triangle **constituents;
	int numberOfConstituentTriangles;

	AbstractTriangle bounder0;
	AbstractTriangle inverse_bounder0;
	AbstractTriangle bounder1;
	AbstractTriangle inverse_bounder1;
	AbstractTriangle bounder2;
	AbstractTriangle inverse_bounder2;
	AbstractTriangle bounder3;
	AbstractTriangle inverse_bounder3;
	AbstractTriangle bounder4;
	AbstractTriangle inverse_bounder4;
	AbstractTriangle bounder5;
	AbstractTriangle inverse_bounder5;
	AbstractTriangle bounder6;
	AbstractTriangle inverse_bounder6;
	AbstractTriangle bounder7;
	AbstractTriangle inverse_bounder7;
	AbstractTriangle bounder8;
	AbstractTriangle inverse_bounder8;
	AbstractTriangle bounder9;
	AbstractTriangle inverse_bounder9;
	AbstractTriangle bounder10;
	AbstractTriangle inverse_bounder10;
	AbstractTriangle bounder11;
	AbstractTriangle inverse_bounder11;
	
	AbstractTriangle *bounders[12];
	AbstractTriangle *inverse_bounders[12];

	float x_min;  
	float x_max;
	float y_min;
	float y_max;
	float z_min;
	float z_max;

	
		
	Triangle** childTriangles[8];

	void create_children(Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike){

		if(numberOfConstituentTriangles <= 300){
			hasChildren = false;
			return;
		}else{
			hasChildren = true;
		}
		children = new Mesh3D*[8];

		float x_seperator = (x_min+x_max)/2;
		float y_seperator = (y_min+y_max)/2;
		float z_seperator = (z_min+z_max)/2;

		int sizes[8] = {0,0,0,0,0,0,0,0};
		int indexes[8] = {0,0,0,0,0,0,0,0};

		for(int i = 0 ; i < numberOfConstituentTriangles; i++){
			Triangle *assesedTriangle = constituents[i];
			
			bool below_x = (assesedTriangle->v1.first <= x_seperator)||
							(assesedTriangle->v2.first <= x_seperator)||
							(assesedTriangle->v3.first <= x_seperator);
			
			bool above_x = (assesedTriangle->v1.first >= x_seperator)||
							(assesedTriangle->v2.first >= x_seperator)||
							(assesedTriangle->v3.first >= x_seperator);
			
			bool below_y = (assesedTriangle->v1.second <= y_seperator)||
							(assesedTriangle->v2.second <= y_seperator)||
							(assesedTriangle->v3.second <= y_seperator);
			
			bool above_y = (assesedTriangle->v1.second >= y_seperator)||
							(assesedTriangle->v2.second >= y_seperator)||
							(assesedTriangle->v3.second >= y_seperator);
			
			bool below_z = (assesedTriangle->v1.third <= z_seperator)||
							(assesedTriangle->v2.third <= z_seperator)||
							(assesedTriangle->v3.third <= z_seperator);
			
			bool above_z = (assesedTriangle->v1.third >= z_seperator)||
							(assesedTriangle->v2.third >= z_seperator)||
							(assesedTriangle->v3.third >= z_seperator);

			if(below_x && below_y && below_z){
				sizes[0]++;	
			}
			
			if(below_x && below_y && above_z){
				sizes[1]++;		
			}

			if(below_x && above_y && below_z){
				sizes[2]++;	
			}
			
			if(below_x && above_y && above_z){
				sizes[3]++;	
			}
			
			if(above_x && below_y && below_z){
				sizes[4]++;	
			}
			
			if(above_x && below_y && above_z){
				sizes[5]++;	
			}
			
			if(above_x && above_y && below_z){
				sizes[6]++;
			}
			
			if(above_x && above_y && above_z){	
				sizes[7]++;
			}

		}

		for(int i = 0 ; i < 8 ; i++){
			if(sizes[i]!=0){
				childTriangles[i] = new Triangle*[sizes[i]];
			}
		}

		
		for(int i = 0 ; i < numberOfConstituentTriangles; i++){
			Triangle *assesedTriangle = constituents[i];
			
			bool below_x = (assesedTriangle->v1.first <= x_seperator)||
							(assesedTriangle->v2.first <= x_seperator)||
							(assesedTriangle->v3.first <= x_seperator);
			
			bool above_x = (assesedTriangle->v1.first >= x_seperator)||
							(assesedTriangle->v2.first >= x_seperator)||
							(assesedTriangle->v3.first >= x_seperator);
			
			bool below_y = (assesedTriangle->v1.second <= y_seperator)||
							(assesedTriangle->v2.second <= y_seperator)||
							(assesedTriangle->v3.second <= y_seperator);
			
			bool above_y = (assesedTriangle->v1.second >= y_seperator)||
							(assesedTriangle->v2.second >= y_seperator)||
							(assesedTriangle->v3.second >= y_seperator);
			
			bool below_z = (assesedTriangle->v1.third <= z_seperator)||
							(assesedTriangle->v2.third <= z_seperator)||
							(assesedTriangle->v3.third <= z_seperator);
			
			bool above_z = (assesedTriangle->v1.third >= z_seperator)||
							(assesedTriangle->v2.third >= z_seperator)||
							(assesedTriangle->v3.third >= z_seperator);

			if(below_x && below_y && below_z){
				childTriangles[0][indexes[0]] = assesedTriangle;
				indexes[0]++;	
			}
			
			if(below_x && below_y && above_z){
				childTriangles[1][indexes[1]] = assesedTriangle;
				indexes[1]++;		
			}

			if(below_x && above_y && below_z){
				childTriangles[2][indexes[2]] = assesedTriangle;
				indexes[2]++;	
			}
			
			if(below_x && above_y && above_z){
				childTriangles[3][indexes[3]] = assesedTriangle;
				indexes[3]++;	
			}
			
			if(above_x && below_y && below_z){
				childTriangles[4][indexes[4]] = assesedTriangle;
				indexes[4]++;	
			}
			
			if(above_x && below_y && above_z){
				childTriangles[5][indexes[5]] = assesedTriangle;
				indexes[5]++;	
			}
			
			if(above_x && above_y && below_z){
				childTriangles[6][indexes[6]] = assesedTriangle;
				indexes[6]++;
			}
			
			if(above_x && above_y && above_z){	
				childTriangles[7][indexes[7]] = assesedTriangle;
				indexes[7]++;
			}

		}

		for(int i = 0 ; i < 8 ; i ++){
			children[i] = new Mesh3D(childTriangles[i],sizes[i], ambientReflectance, 
			 diffuseReflectance,
			 specularReflectance, 
			 mirrorReflectance ,
			 phongExp,
			 isMirrorLike);
			
		}


	}

	

	void createBounders(float peculiar_a0, float peculiar_a1, float peculiar_b0, float peculiar_b1 ,float common, int common_index,int pairIndex){
		Vector3f v0;
		Vector3f v1;
		Vector3f v2;
		Vector3f v3;

		if(common_index == 0){
			v0.second = peculiar_a0; 
			v0.third = peculiar_b0;
			
			v1.second = peculiar_a0;
			v1.third = peculiar_b1;

			v2.second = peculiar_a1;
			v2.third = peculiar_b1;
			
			v3.second = peculiar_a1;
			v3.third = peculiar_b0;

			v0.first = common;
			v1.first = common;
			v2.first = common;
			v3.first = common;


		}else if(common_index == 1){
			v0.first = peculiar_a0; 
			v0.third = peculiar_b0;
			
			v1.first = peculiar_a0;
			v1.third = peculiar_b1;

			v2.first = peculiar_a1;
			v2.third = peculiar_b1;
			
			v3.first = peculiar_a1;
			v3.third = peculiar_b0;


			v0.second = common;
			v1.second = common;
			v2.second = common;
			v3.second = common;


		}else if(common_index == 2){
			v0.first = peculiar_a0; 
			v0.second = peculiar_b0;
			
			v1.first = peculiar_a0;
			v1.second = peculiar_b1;

			v2.first = peculiar_a1;
			v2.second = peculiar_b1;
			
			v3.first = peculiar_a1;
			v3.second = peculiar_b0;


			v0.third = common;
			v1.third = common;
			v2.third = common;
			v3.third = common;
		}
		
		bounders[2*pairIndex]->setup(v0,v1,v3); 
		bounders[2*pairIndex+1]->setup(v1,v2,v3); 
		
		
		inverse_bounders[2*pairIndex]->setup(v3,v1,v0); 
		inverse_bounders[2*pairIndex+1]->setup(v3,v2,v1); 
		

	}

	public:
	Mesh3D(parser::Face *faces, int numOfFaces, 	
			Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike
	):bounder0(),bounder1(),bounder2(),bounder3(),bounder4(),bounder5(),bounder6(),bounder7(),bounder8(),bounder9(),bounder10(),bounder11(),
	inverse_bounder0(),inverse_bounder1(),inverse_bounder2(),inverse_bounder3(),inverse_bounder4(),inverse_bounder5(),inverse_bounder6(),inverse_bounder7(),inverse_bounder8(),inverse_bounder9(),inverse_bounder10(),inverse_bounder11(),
			Object3D(
					ambientReflectance, 
					diffuseReflectance,
					specularReflectance, 
					mirrorReflectance ,
					phongExp,
					isMirrorLike
				){
					bounders[0]= &bounder0;
					inverse_bounders[0]= &inverse_bounder0;
					bounders[1]= &bounder1;
					inverse_bounders[1]= &inverse_bounder1;
					bounders[2]= &bounder2;
					inverse_bounders[2]= &inverse_bounder2;
					bounders[3]= &bounder3;
					inverse_bounders[3]= &inverse_bounder3;
					bounders[4]= &bounder4;
					inverse_bounders[4]= &inverse_bounder4;
					bounders[5]= &bounder5;
					inverse_bounders[5]= &inverse_bounder5;
					bounders[6]= &bounder6;
					inverse_bounders[6]= &inverse_bounder6;
					bounders[7]= &bounder7;
					inverse_bounders[7]= &inverse_bounder7;
					bounders[8]= &bounder8;
					inverse_bounders[8]= &inverse_bounder8;
					bounders[9]= &bounder9;
					inverse_bounders[9]= &inverse_bounder9;
					bounders[10]= &bounder10;
					inverse_bounders[10]= &inverse_bounder10;
					bounders[11]= &bounder11;
					inverse_bounders[11]= &inverse_bounder11;
					
					constituents = new Triangle*[numOfFaces];
					numberOfConstituentTriangles = numOfFaces;
					bool instanciated = false;

					for(int i = 0 ; i < numOfFaces;i++){
						Vector3f first_vertex = getVertex(faces[i].v0_id);
						Vector3f second_vertex = getVertex(faces[i].v1_id);
						Vector3f third_vertex = getVertex(faces[i].v2_id);
						
						if(!instanciated){
							x_min = min(first_vertex.first ,min(second_vertex.first,third_vertex.first));
							x_max = max(first_vertex.first ,max(second_vertex.first,third_vertex.first));

							y_min = min(first_vertex.second,min(second_vertex.second,third_vertex.second));
							y_max = max(first_vertex.second,max(second_vertex.second,third_vertex.second));

							z_min = min(first_vertex.third,min(second_vertex.third,third_vertex.third));
							z_max = max(first_vertex.third,max(second_vertex.third,third_vertex.third));
							instanciated = true;
						}
						
						x_min = min(first_vertex.first, x_min);			
						x_max = max(first_vertex.first, x_max);
						x_min = min(second_vertex.first, x_min);			
						x_max = max(second_vertex.first, x_max);	
						x_min = min(third_vertex.first, x_min);			
						x_max = max(third_vertex.first, x_max);
						
						y_min = min(first_vertex.second, y_min);			
						y_max = max(first_vertex.second, y_max);
						y_min = min(second_vertex.second, y_min);			
						y_max = max(second_vertex.second, y_max);	
						y_min = min(third_vertex.second, y_min);			
						y_max = max(third_vertex.second, y_max);
						
						z_min = min(first_vertex.third, z_min);			
						z_max = max(first_vertex.third, z_max);
						z_min = min(second_vertex.third, z_min);			
						z_max = max(second_vertex.third, z_max);	
						z_min = min(third_vertex.third, z_min);			
						z_max = max(third_vertex.third, z_max);							
					


						constituents[i] = new Triangle(first_vertex,
														second_vertex,
														third_vertex,
														ambientReflectance,
														diffuseReflectance,
														specularReflectance,
														mirrorReflectance,
														phongExp,
														isMirrorLike);
					}
				
				if(numberOfConstituentTriangles != 0){
						
					createBounders(x_min,x_max,y_min,y_max,z_max,2,0);
					createBounders(x_min,x_max,y_min,y_max,z_min,2,1);
					
					createBounders(x_min,x_max,z_min,z_max,y_max,1,2);
					createBounders(x_min,x_max,z_min,z_max,y_min,1,3);
					
					createBounders(y_min,y_max,z_min,z_max,x_max,0,4);
					createBounders(y_min,y_max,z_min,z_max,x_min,0,5);
					create_children(ambientReflectance, 
									diffuseReflectance,
									specularReflectance, 
									mirrorReflectance ,
									phongExp,
									isMirrorLike);
				}
	}
	Mesh3D(Triangle **triangles, int numberOfTriangles, 	
			Vector3f ambientReflectance, 
			Vector3f diffuseReflectance,
			Vector3f specularReflectance, 
			Vector3f mirrorReflectance ,
			int phongExp,
			bool isMirrorLike
	):bounder0(),bounder1(),bounder2(),bounder3(),bounder4(),bounder5(),bounder6(),bounder7(),bounder8(),bounder9(),bounder10(),bounder11(),
	inverse_bounder0(),inverse_bounder1(),inverse_bounder2(),inverse_bounder3(),inverse_bounder4(),inverse_bounder5(),inverse_bounder6(),inverse_bounder7(),inverse_bounder8(),inverse_bounder9(),inverse_bounder10(),inverse_bounder11(),
	 		Object3D(
					ambientReflectance, 
					diffuseReflectance,
					specularReflectance, 
					mirrorReflectance ,
					phongExp,
					isMirrorLike
				){
					
					bounders[0]= &bounder0;
					inverse_bounders[0]= &inverse_bounder0;
					
					bounders[1]= &bounder1;
					inverse_bounders[1]= &inverse_bounder1;
					
					bounders[2]= &bounder2;
					inverse_bounders[2]= &inverse_bounder2;
					
					bounders[3]= &bounder3;
					inverse_bounders[3]= &inverse_bounder3;
					
					bounders[4]= &bounder4;
					inverse_bounders[4]= &inverse_bounder4;
					
					bounders[5]= &bounder5;
					inverse_bounders[5]= &inverse_bounder5;
					
					bounders[6]= &bounder6;
					inverse_bounders[6]= &inverse_bounder6;
					
					bounders[7]= &bounder7;
					inverse_bounders[7]= &inverse_bounder7;
					
					bounders[8]= &bounder8;
					inverse_bounders[8]= &inverse_bounder8;
					
					bounders[9]= &bounder9;
					inverse_bounders[9]= &inverse_bounder9;
					
					bounders[10]= &bounder10;
					inverse_bounders[10]= &inverse_bounder10;
					
					bounders[11]= &bounder11;
					inverse_bounders[11]= &inverse_bounder11;
					

					
					constituents = triangles;
					numberOfConstituentTriangles = numberOfTriangles;
					bool instanciated = false;

					for(int i = 0 ; i < numberOfTriangles;i++){
						Vector3f first_vertex = triangles[i]->v1 ;
						Vector3f second_vertex = triangles[i]->v2;
						Vector3f third_vertex = triangles[i]->v3;
						
						if(!instanciated){
							x_min = min(first_vertex.first ,min(second_vertex.first,third_vertex.first));
							x_max = max(first_vertex.first ,max(second_vertex.first,third_vertex.first));

							y_min = min(first_vertex.second,min(second_vertex.second,third_vertex.second));
							y_max = max(first_vertex.second,max(second_vertex.second,third_vertex.second));

							z_min = min(first_vertex.third,min(second_vertex.third,third_vertex.third));
							z_max = max(first_vertex.third,max(second_vertex.third,third_vertex.third));
							instanciated = true;
						}
						
						x_min = min(first_vertex.first, x_min);			
						x_max = max(first_vertex.first, x_max);
						x_min = min(second_vertex.first, x_min);			
						x_max = max(second_vertex.first, x_max);	
						x_min = min(third_vertex.first, x_min);			
						x_max = max(third_vertex.first, x_max);
						
						y_min = min(first_vertex.second, y_min);			
						y_max = max(first_vertex.second, y_max);
						y_min = min(second_vertex.second, y_min);			
						y_max = max(second_vertex.second, y_max);	
						y_min = min(third_vertex.second, y_min);			
						y_max = max(third_vertex.second, y_max);
						
						z_min = min(first_vertex.third, z_min);			
						z_max = max(first_vertex.third, z_max);
						z_min = min(second_vertex.third, z_min);			
						z_max = max(second_vertex.third, z_max);	
						z_min = min(third_vertex.third, z_min);			
						z_max = max(third_vertex.third, z_max);							
					}
				
				if(numberOfConstituentTriangles != 0){

					createBounders(x_min,x_max,y_min,y_max,z_max,2,0);
					createBounders(x_min,x_max,y_min,y_max,z_min,2,1);
					
					createBounders(x_min,x_max,z_min,z_max,y_max,1,2);
					createBounders(x_min,x_max,z_min,z_max,y_min,1,3);
					
					createBounders(y_min,y_max,z_min,z_max,x_max,0,4);
					createBounders(y_min,y_max,z_min,z_max,x_min,0,5);
					create_children(ambientReflectance, 
									diffuseReflectance,
									specularReflectance, 
									mirrorReflectance ,
									phongExp,
									isMirrorLike);


				}
				
		}
	
	
	
	
	
	
	
	
	
	virtual IntersectionReport getIntersectionReport(Ray ray) override{
			
		IntersectionReport result;
		if (numberOfConstituentTriangles == 0){
			result.occurance = false;
			return result;
		}
		
		bool include = false;
		
		for(int i = 0 ; i < 12 ; i++){
			if(bounders[i]->getIntersectionReport(ray).occurance || inverse_bounders[i]->getIntersectionReport(ray).occurance ){
				include = true;
				break;
			}
		}
		
		if(!include){
			result.occurance = false;
			return result;
		}
		if(hasChildren){
			result.occurance = false;

			for(int i = 0 ; i < 8 ; i++){
				IntersectionReport temporaryResult = children[i]->getIntersectionReport(ray);
				
				if(result.occurance){
					if(temporaryResult.occurance && temporaryResult.t_parameter < result.t_parameter){
						result = temporaryResult;
					}
				}else{
					result = temporaryResult;
				}
			}

			return result;
		}

		include = false;
		IntersectionReport effectiveIntersection;
		float min_t_parameter ; 

		for(int i= 0 ; i < numberOfConstituentTriangles; i ++ ){
			IntersectionReport acquiredReport = constituents[i]->getIntersectionReport(ray);
			if(acquiredReport.occurance){
				if(!include){
					effectiveIntersection = acquiredReport;
					include = true;
					min_t_parameter = acquiredReport.t_parameter;
				}else if(acquiredReport.t_parameter < min_t_parameter){
					effectiveIntersection = acquiredReport;
					min_t_parameter = acquiredReport.t_parameter;
				}
			}
		}


		if(include){
			return effectiveIntersection;
		}else{
			result.occurance = false;
			return result;
		}
		
	}
};


/*globals*/





Vector3f computeRayColor(Ray ray, int limit){
	if(limit == 0){
		throw "ERROR";
	}

	int objectIndex;
	bool objectHit = false;
	float min_t_parameter;
	IntersectionReport acquiredIntersectionReport;
	for(int i = 0 ; i < object_count ; i++){
		IntersectionReport intersectionReport = objects[i]->getIntersectionReport(ray);
		if(intersectionReport.occurance){
				if(!objectHit || intersectionReport.t_parameter < min_t_parameter){
					min_t_parameter = intersectionReport.t_parameter;
					acquiredIntersectionReport = intersectionReport;
					objectHit = true;
					objectIndex = i;
				}
		}
	}

	Vector3f result;
	if(objectHit){
		result = objects[objectIndex]->getVisual(acquiredIntersectionReport,limit -1);
	}else{
		result.first = scene.background_color.x;
		result.second = scene.background_color.y;
		result.third = scene.background_color.z;
	}
	return result;
}

class Camera{
    private:
	Vector3f pixel_0_0_coordinate;
	Vector3f lookRight;
	Vector3f gaze;
	Vector3f lookUp;
	Vector3f viewPoint; 
	float pixel_width;
	float pixel_height;
	
	Vector3f computePixelCoordinate(Vector2i pixel){
		Vector3f result = this->pixel_0_0_coordinate;
		result = evaluator.add(result,evaluator.multiply(this->lookRight,pixel.first * this->pixel_width));
		result = evaluator.substract(result,evaluator.multiply(this->lookUp,pixel.second * this->pixel_height));
		
		return result;
	}
	
	

    Vector3i computePixelColor(unsigned x, unsigned y){

		Ray ray;
		ray.source = viewPoint;
		Vector2i pixel;
		pixel.first = x;
		pixel.second = y;
		ray.direction= evaluator.normalize(evaluator.substract(computePixelCoordinate(pixel),this->viewPoint));

		Vector3f result = computeRayColor(ray,scene.max_recursion_depth+1);
		return evaluator.toIntVectorByLimit(result,255);
    }
    parser::Camera descriptor;

    public:
    Camera(parser::Camera descriptor){
        this->descriptor = descriptor;
		gaze.first = descriptor.gaze.x;
		gaze.second = descriptor.gaze.y;
		gaze.third = descriptor.gaze.z;
		lookUp.first = descriptor.up.x;
		lookUp.second = descriptor.up.y;
		lookUp.third = descriptor.up.z;

		this->lookRight = evaluator.crossProduct(this->gaze,this->lookUp);
		this->pixel_width = (descriptor.near_plane.y - descriptor.near_plane.x)/descriptor.image_width;
		this->pixel_height = (descriptor.near_plane.w - descriptor.near_plane.z)/descriptor.image_height;
		
		
		this->viewPoint.first = descriptor.position.x;
		this->viewPoint.second = descriptor.position.y;
		this->viewPoint.third = descriptor.position.z;

		
		Vector3f center = evaluator.add(this->viewPoint,evaluator.multiply(this->gaze,descriptor.near_distance));
		Vector3f topLeft = evaluator.add(evaluator.add(center,evaluator.multiply(this->lookUp,descriptor.near_plane.w)),evaluator.multiply(this->lookRight,descriptor.near_plane.x));	
		
		Vector3f temp = topLeft;
		temp.first += this->pixel_width/2;
		temp.second -= this->pixel_height/2;	
		this->pixel_0_0_coordinate = temp;
    }
	


    void computeImage(){
        unsigned char* image = new unsigned char [descriptor.image_width * descriptor.image_height * 3];

        int i = 0;
        for (int y = 0; y < descriptor.image_height; ++y)
        {
            for (int x = 0; x < descriptor.image_width; ++x)
            {
				
                Vector3i pixelColor = computePixelColor(x,y);
                image[i++] = pixelColor.first;
                image[i++] = pixelColor.second;
                image[i++] = pixelColor.third;
            }
        }
        write_ppm(descriptor.image_name.c_str(), image, descriptor.image_width, descriptor.image_height);
    }
};

struct Material{
	bool isMirrorLike;
	Vector3f ambientReflectance; 
	Vector3f diffuseReflectance;
	Vector3f specularReflectance; 
	Vector3f mirrorReflectance;
	int phongExp;
			
};



int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    
	
    scene.loadFromXml(argv[1]);

	//populate pointLightSources,pointLigtSurceCount
	
	pointLightSourceCount = scene.point_lights.size();
	pointLightSources = new PointLightSource[pointLightSourceCount];
	for(int i=0; i<pointLightSourceCount; i++){
		pointLightSources[i].position = fromParserToUsVector3fConverter(scene.point_lights[i].position);
		pointLightSources[i].color = fromParserToUsVector3fConverter(scene.point_lights[i].intensity);
	}

	
	int sphereSize = scene.spheres.size();
	int triangleSize = scene.triangles.size();
	int meshSize = scene.meshes.size();

	object_count = sphereSize + triangleSize + meshSize;
	objects = new Object3D*[object_count];

	int i = 0 ;
	for(int j = 0 ; j < sphereSize; j++){
		parser::Sphere currentSphereFromParser = scene.spheres[j];
		parser::Vec3f currentCenterFromParser = scene.vertex_data[currentSphereFromParser.center_vertex_id-1];
		float currentRadiusFromParser = currentSphereFromParser.radius;

		parser::Material currentMaterialFromParser = scene.materials[currentSphereFromParser.material_id-1];

		Sphere *currentSphere = new Sphere(fromParserToUsVector3fConverter(currentCenterFromParser),
										currentRadiusFromParser,
										fromParserToUsVector3fConverter(currentMaterialFromParser.ambient),
										fromParserToUsVector3fConverter(currentMaterialFromParser.diffuse),
										fromParserToUsVector3fConverter(currentMaterialFromParser.specular),
										fromParserToUsVector3fConverter(currentMaterialFromParser.mirror),
										currentMaterialFromParser.phong_exponent,
										currentMaterialFromParser.is_mirror
										);
		objects[i] = (Object3D *) currentSphere;
		i +=1;
	}

	for(int j = 0 ; j < triangleSize; j++){
		parser::Triangle currentTriangleFromParser = scene.triangles[j];
		parser::Vec3f indice0 = scene.vertex_data[currentTriangleFromParser.indices.v0_id-1];
		parser::Vec3f indice1 = scene.vertex_data[currentTriangleFromParser.indices.v1_id-1];
		parser::Vec3f indice2 = scene.vertex_data[currentTriangleFromParser.indices.v2_id-1];

		parser::Material currentMaterialFromParser = scene.materials[currentTriangleFromParser.material_id-1];

		Triangle *currentTriangle = new Triangle(fromParserToUsVector3fConverter(indice0),
											fromParserToUsVector3fConverter(indice1),
											fromParserToUsVector3fConverter(indice2),
												fromParserToUsVector3fConverter(currentMaterialFromParser.ambient),
												fromParserToUsVector3fConverter(currentMaterialFromParser.diffuse),
												fromParserToUsVector3fConverter(currentMaterialFromParser.specular),
												fromParserToUsVector3fConverter(currentMaterialFromParser.mirror),
												currentMaterialFromParser.phong_exponent,
												currentMaterialFromParser.is_mirror
												);
		objects[i] = (Object3D *) currentTriangle;
		i +=1;
	}
	for(int j = 0 ; j < meshSize; j++){
		parser::Mesh currentMeshFromParser = scene.meshes[j];

		parser::Material currentMaterialFromParser = scene.materials[currentMeshFromParser.material_id-1];

		Mesh3D *currentMesh =new Mesh3D(currentMeshFromParser.faces.data(),currentMeshFromParser.faces.size() ,
										fromParserToUsVector3fConverter(currentMaterialFromParser.ambient),
										fromParserToUsVector3fConverter(currentMaterialFromParser.diffuse),
										fromParserToUsVector3fConverter(currentMaterialFromParser.specular),
										fromParserToUsVector3fConverter(currentMaterialFromParser.mirror),
										currentMaterialFromParser.phong_exponent,
										currentMaterialFromParser.is_mirror
										);
		objects[i] = (Object3D *)  currentMesh;
		i +=1;
	}

	


    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    for(int i = 0 ; i < scene.cameras.size(); i++){
        Camera camera(scene.cameras[i]);
        camera.computeImage();
    }

	
	return 0;
  
   

    
}