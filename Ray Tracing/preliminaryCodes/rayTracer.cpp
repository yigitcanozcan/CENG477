#include<iostream>
#include<math.h>
using namespace std;



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


class LightSource{
	public:
	virtual int getLightSourceType(){
		return -1;
	}
};
class PointLightSource : LightSource{
	public:
	Vector3f color ;
	Vector3f position;
	PointLightSource(Vector3f position_parameter, Vector3f color_parameter){
		color = color_parameter;
		position = position_parameter;
	}
	virtual int getLightSourceType() override{
		return 0;
	}
};
Vector3f computeRayColorWithObjectProhibition(Ray ray, int limit, int prohibitedObjectIndex);

class Object3D;
struct {
	struct
	{
		int width;
		int height;
	}resolution;
	Vector3i bgColor;
	struct
	{
		Vector3f viewPoint;
		Vector3f gaze;
		Vector3f lookUp;
		int focalLength;
		struct 
		{
			float l,r,t,b;
		} viewWindow;
		
	} viewStruct;
	struct{
		int objectCount;
		Object3D **objects;
	} objects_3d;
	struct
	{
		int lightSourceCount;
		LightSource **sources;
	} lightSources;
	Vector3f ambientLight;

} confStruct;
class Object3D{
	public:
	virtual IntersectionReport getIntersectionReport(Ray ray){throw;}
	int index;
	float reflectanceAttenuationCoefficient;
	int specularCoefficient;
	Vector3f color;
	Object3D(Vector3f color_parameter, int specularCoefficient_parameter, 
				float reflectanceAttenuationCoefficient_parameter, int index_parameter){
			color = color_parameter;
			specularCoefficient = specularCoefficient_parameter;
			reflectanceAttenuationCoefficient = reflectanceAttenuationCoefficient_parameter;
			index = index_parameter;
	}
	virtual Vector3f getVisual(IntersectionReport intersectionReport, int limit){	
		//ambient component
		Vector3f result = evaluator.elementwiseMultiplication(confStruct.ambientLight,color);
		for(int i = 0 ; i < confStruct.lightSources.lightSourceCount ; i++){
			LightSource *currentSource = confStruct.lightSources.sources[i];
			if(currentSource->getLightSourceType()==0){
				PointLightSource *currentPointLS = (PointLightSource *) currentSource;
				
				//shadowing
				bool prematureIntersection = false;
				for(int j = 0 ; j < confStruct.objects_3d.objectCount; j++){
					if(j == this->index){
						continue;
					}
					Ray computedRay;
					computedRay.source = intersectionReport.intersectionPoint;
					computedRay.direction = evaluator.substract(currentPointLS->position, intersectionReport.intersectionPoint);
					IntersectionReport computedIntersection = confStruct.objects_3d.objects[j]->getIntersectionReport(computedRay);
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

				
				
				
				
				
				Vector3f illuminationVector = evaluator.normalize(evaluator.substract(currentPointLS->position,intersectionReport.intersectionPoint));
				//diffuse component
				float diffusionCoefficient = evaluator.dotProduct(illuminationVector,intersectionReport.intersectionSurfaceNormal);
				if(diffusionCoefficient > 0){
					Vector3f contribution = evaluator.elementwiseMultiplication(evaluator.multiply(currentPointLS->color,diffusionCoefficient),color);
					result = evaluator.add(result,contribution);
				}
				//specular component
				Vector3f to_observer = evaluator.normalize(evaluator.substract(intersectionReport.intersectedRay.source,intersectionReport.intersectionPoint));
				Vector3f halfVector = evaluator.normalize(evaluator.add(to_observer,illuminationVector));
				float temp = evaluator.dotProduct(halfVector, intersectionReport.intersectionSurfaceNormal);
				if(temp > 0){
					float specularContributionCoefficient = pow(temp, specularCoefficient);
					Vector3f specularContribution = evaluator.multiply(currentPointLS->color, specularContributionCoefficient);
					result = evaluator.add(result,specularContribution);
				}

			}else{
				throw;
			}
		}
		//computation of reflectance
		if((limit > 0) && (this->reflectanceAttenuationCoefficient > 0)){
			Ray computationRequest;
			computationRequest.source = intersectionReport.intersectionPoint;
			computationRequest.direction = evaluator.add(intersectionReport.intersectedRay.direction, 
				evaluator.multiply(intersectionReport.intersectionSurfaceNormal,
					(-2)*evaluator.dotProduct(intersectionReport.intersectedRay.direction,intersectionReport.intersectionSurfaceNormal)
				)
			);
			Vector3f computedRayColor = computeRayColorWithObjectProhibition(computationRequest, limit-1, this->index);
			Vector3f computedRayColorContribution = evaluator.multiply(computedRayColor, this->reflectanceAttenuationCoefficient);
			result = evaluator.add(result, computedRayColorContribution);
		}
		return result;
		
	}
};

Vector3f computeRayColorWithObjectProhibition(Ray ray, int limit, int prohibitedObjectIndex){
	//TODO: Ambiance or bgcolor
	Vector3f result = evaluator.toFloatVector(confStruct.bgColor) ;
	if(limit == 0){
		return result;
	}
	int objectIndex;
	bool objectHit = false;
	float min_t_parameter;
	IntersectionReport acquiredIntersectionReport;
	for(int i = 0 ; i < confStruct.objects_3d.objectCount ; i++){
		if(i == prohibitedObjectIndex){
			continue;
		}
		IntersectionReport intersectionReport = confStruct.objects_3d.objects[i]->getIntersectionReport(ray);
		if(intersectionReport.occurance){
				if(!objectHit || intersectionReport.t_parameter < min_t_parameter){
					min_t_parameter = intersectionReport.t_parameter;
					acquiredIntersectionReport = intersectionReport;
					objectHit = true;
					objectIndex = i;
				}
		}
	}

	
	if(objectHit){
		result = confStruct.objects_3d.objects[objectIndex]->getVisual(acquiredIntersectionReport,limit -1 );
	}
	return result;
}

class Sphere :public Object3D{
	public:
	float radius;
	Vector3f center;
	
	Sphere(Vector3f parameter_center, float parameter_radius, Vector3f parameter_color, int specularCoeficientParameter,
			float reflectanceAttenuationCoefficient_parameter, int index_parameter): 
				Object3D(parameter_color,specularCoeficientParameter,
					reflectanceAttenuationCoefficient_parameter, index_parameter){
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

class Triangle: public Object3D{
	public:
	Vector3f v1,v2,v3;
	Triangle(Vector3f v1_parameter, Vector3f v2_parameter, Vector3f v3_parameter, Vector3f parameter_color, int specularCoeficientParameter,
			float reflectanceAttenuationCoefficient_parameter, int index_parameter): 
				Object3D(parameter_color,specularCoeficientParameter,
					reflectanceAttenuationCoefficient_parameter, index_parameter){
						v1 = v1_parameter;
						v2 = v2_parameter;
						v3 = v3_parameter;
	}
	virtual IntersectionReport getIntersectionReport(Ray ray) override{
		IntersectionReport result;
		result.intersectedRay = ray;
		result.intersectionSurfaceNormal = evaluator.normalize(
			evaluator.crossProduct(
				evaluator.substract(v2,v1),
				evaluator.substract(v3,v1)
			)
		);
		result.occurance = true;
		float c1,c2,t;
		Matrix3x3 matrix = {
			evaluator.substract(v1,v3),
			evaluator.substract(v2,v3),
			evaluator.multiply(ray.direction,-1)
		};
		float determinant = evaluator.computeDeterminant3x3(matrix);
		
		if(determinant == 0){
			result.occurance = false;
		}else{
			Vector3f resultVector = evaluator.substract(ray.source,v3);
			
			Matrix3x3 activeMatrixC1 = matrix;
			activeMatrixC1.col0 = resultVector;
			c1 = evaluator.computeDeterminant3x3(activeMatrixC1)/determinant;

			Matrix3x3 activeMatrixC2 = matrix;
			activeMatrixC2.col1 = resultVector;
			c2 = evaluator.computeDeterminant3x3(activeMatrixC2)/determinant;

			Matrix3x3 activeMatrixT = matrix;
			activeMatrixT.col2 = resultVector;
			t = evaluator.computeDeterminant3x3(activeMatrixT)/determinant;
			if(t < 0){
				result.occurance = false;
			}
			if(c1 < 0 || c1 > 1){
				result.occurance = false;
			}
			if(c2 < 0 || c2 > 1){
				result.occurance = false;
			}
			if(1-c1-c2 < 0 || 1-c1-c2 > 1){
				result.occurance = false;
			}
		}
		
		if(result.occurance){
			result.t_parameter = t; 
			result.intersectionPoint = evaluator.add(ray.source,
				evaluator.multiply(ray.direction,t)
			);
		}
		return result;
	}
};





Vector3f computeRayColor(Ray ray, int limit = 3){
	Vector3f result = evaluator.toFloatVector(confStruct.bgColor) ;
	if(limit == 0){
		return result;
	}
	int objectIndex;
	bool objectHit = false;
	float min_t_parameter;
	IntersectionReport acquiredIntersectionReport;
	for(int i = 0 ; i < confStruct.objects_3d.objectCount ; i++){
		IntersectionReport intersectionReport = confStruct.objects_3d.objects[i]->getIntersectionReport(ray);
		if(intersectionReport.occurance){
				if(!objectHit || intersectionReport.t_parameter < min_t_parameter){
					min_t_parameter = intersectionReport.t_parameter;
					acquiredIntersectionReport = intersectionReport;
					objectHit = true;
					objectIndex = i;
				}
		}
	}

	
	if(objectHit){
		result = confStruct.objects_3d.objects[objectIndex]->getVisual(acquiredIntersectionReport,limit -1 );
	}
	return result;
}

struct {
	Vector3f lookRight;
	Vector3f pixel_0_0_coordinate;
	float pixel_width;
	float pixel_height;
} helperStruct ;

Vector3f computePixelCoordinate(Vector2i pixel){
	Vector3f result = helperStruct.pixel_0_0_coordinate;
	result = evaluator.add(result,evaluator.multiply(helperStruct.lookRight,pixel.second*helperStruct.pixel_width));
	result = evaluator.substract(result,evaluator.multiply(confStruct.viewStruct.lookUp,pixel.first*helperStruct.pixel_height));
	
	return result;
}

Vector3i computePixelColor(Vector2i pixel){
	Ray ray;
	Vector3f pixelCoordinate = computePixelCoordinate(pixel);
	ray.direction = evaluator.substract(pixelCoordinate,confStruct.viewStruct.viewPoint);
	
	ray.source = confStruct.viewStruct.viewPoint;
	return evaluator.toIntVectorByLimit(evaluator.multiply(computeRayColor(ray),255) ,255);
}
void populateHelperStruct(){
	helperStruct.lookRight = evaluator.crossProduct(confStruct.viewStruct.gaze,confStruct.viewStruct.lookUp);
	helperStruct.pixel_width = (confStruct.viewStruct.viewWindow.r-confStruct.viewStruct.viewWindow.l)/confStruct.resolution.width;
	helperStruct.pixel_height = (confStruct.viewStruct.viewWindow.t-confStruct.viewStruct.viewWindow.b)/confStruct.resolution.height;
	Vector3f center = evaluator.add(confStruct.viewStruct.viewPoint,evaluator.multiply(confStruct.viewStruct.gaze,confStruct.viewStruct.focalLength));
	Vector3f topLeft = evaluator.add(evaluator.add(center,evaluator.multiply(confStruct.viewStruct.lookUp,confStruct.viewStruct.viewWindow.t)),evaluator.multiply(helperStruct.lookRight,confStruct.viewStruct.viewWindow.l));	
	
	Vector3f temp = topLeft;
	temp.first += helperStruct.pixel_width/2;
	temp.second -= helperStruct.pixel_height/2;	
	helperStruct.pixel_0_0_coordinate = temp;
}
void parse(){
	cin >> confStruct.resolution.width >> confStruct.resolution.height  ;
	cin >> confStruct.bgColor.first >> confStruct.bgColor.second >> confStruct.bgColor.third;
	cin >> confStruct.viewStruct.viewPoint.first >> confStruct.viewStruct.viewPoint.second >> confStruct.viewStruct.viewPoint.third;
	cin >> confStruct.viewStruct.gaze.first >> confStruct.viewStruct.gaze.second >> confStruct.viewStruct.gaze.third; 
	cin >> confStruct.viewStruct.lookUp.first >> confStruct.viewStruct.lookUp.second >> confStruct.viewStruct.lookUp.third;
	
	
	confStruct.viewStruct.lookUp = evaluator.normalize(confStruct.viewStruct.lookUp);
	confStruct.viewStruct.gaze = evaluator.normalize(confStruct.viewStruct.gaze);

	cin >> confStruct.viewStruct.viewWindow.l >> confStruct.viewStruct.viewWindow.r >> confStruct.viewStruct.viewWindow.t >> confStruct.viewStruct.viewWindow.b ;
	cin >> confStruct.viewStruct.focalLength;
	cin >> confStruct.objects_3d.objectCount;
	
	if(confStruct.objects_3d.objectCount != 0){
		confStruct.objects_3d.objects = new Object3D*[confStruct.objects_3d.objectCount];
	}
	for(int i = 0 ; i < confStruct.objects_3d.objectCount ; i++){
		int objectType ;
		int index;
		cin >> objectType;
		if (objectType==0){
			float x,y,z,r;
			float red,green,blue;
			int specularCoefficient;
			float reflectanceAttenuationCoefficient;

			cin >> x >> y >> z >> r;
			cin >> specularCoefficient;
			cin >> reflectanceAttenuationCoefficient;
			cin >> red >> green >> blue;

			Vector3f center_position = {x,y,z};
			Vector3f color = {red,green,blue};

			Sphere *createdSphere = new Sphere(center_position, r , color, specularCoefficient,reflectanceAttenuationCoefficient,i);
			confStruct.objects_3d.objects[i] = (Object3D*)createdSphere;
		}else if(objectType==1){
			Vector3f v1,v2,v3;
			float red,green,blue;
			int specularCoefficient;
			float reflectanceAttenuationCoefficient;

			cin >>  v1.first >> v1.second >> v1.third >> 
					v2.first >> v2.second >> v2.third >>
					v3.first >> v3.second >> v3.third ;
			cin >> specularCoefficient;
			cin >> reflectanceAttenuationCoefficient;
			cin >> red >> green >> blue;

			Vector3f color = {red,green,blue};

			Triangle *createdTriangle = new Triangle(v1,v2,v3 , color, specularCoefficient,reflectanceAttenuationCoefficient,i);
			confStruct.objects_3d.objects[i] = (Object3D*)createdTriangle;

		}else{
			throw;
		}
	}
	cin >> confStruct.ambientLight.first >> confStruct.ambientLight.second >> confStruct.ambientLight.third;

	cin >> confStruct.lightSources.lightSourceCount;
	if(confStruct.lightSources.lightSourceCount > 0){
		confStruct.lightSources.sources = new LightSource*[confStruct.lightSources.lightSourceCount];
	}
	for(int i = 0 ; i < confStruct.lightSources.lightSourceCount ; i++){
		int lightSourceType ;
		cin >> lightSourceType;
		if(lightSourceType == 0){
			float x,y,z,r,g,b;
			cin >> x >> y >> z >> r >> g >> b ; 
			confStruct.lightSources.sources[i] = (LightSource *) new PointLightSource({x,y,z},{r,g,b});
		}else{
			throw;
		}
	}
}
void writePpm(Vector3i **image){
	int width = confStruct.resolution.width;
	int height = confStruct.resolution.height; 
	
	cout << "P3" << endl;
	cout << width << " " << height << endl;
	cout << "255" << endl ; 

	for(int i = 0 ; i < height ; i++){
		for(int j = 0 ; j < width; j++){
			cout << image[i][j].first;
			cout << " ";
			cout << image[i][j].second;
			cout << " ";
			cout << image[i][j].third;
			cout << " ";
		}
		cout << endl;
	}
	cout << (unsigned char) 15;
	cout << (unsigned char) 0 ;
	cout << (unsigned char) 0 ; 


	cout << (unsigned char) 0 ;
	cout << (unsigned char) 15;
	cout << (unsigned char) 0 ;
	
	cout << (unsigned char) 0 ; 
	cout << (unsigned char) 0 ;
	cout << (unsigned char) 15; 
}
int main(){
	parse();
	populateHelperStruct();
	int height = confStruct.resolution.height;
	int width = confStruct.resolution.width;

	
	Vector3i **image = new Vector3i* [height];
	for(int i = 0 ; i < height ; i++){
		image[i] = new Vector3i [width];
	}
	
	for(int i = 0 ; i < height ; i++){
		for(int j = 0 ; j < width; j++){
			Vector2i pixel = {i,j};
			image[i][j] = computePixelColor(pixel);
		}
	}

	writePpm(image);

}
