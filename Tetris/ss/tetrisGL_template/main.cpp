//
// Author: Ahmet Oguz Akyuz
// 
// This is a sample code that draws a single block piece at the center
// of the window. It does many boilerplate work for you -- but no
// guarantees are given about the optimality and bug-freeness of the
// code. You can use it to get started or you can simply ignore its
// presence. You can modify it in any way that you like.
//
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <ft2build.h>
#include FT_FREETYPE_H
#include <random>

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

//behaviour configuration
int turningGraduality = 20;
string sideNames[4] = {"front","right","back","left"};
int sideIndex = 0;




GLuint gProgram[4];
int gWidth = 600, gHeight = 1000;
GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLuint gTex2D;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
int gTriangleIndexDataSizeInBytes, gLineIndexDataSizeInBytes;

unsigned state[9][9];
unsigned proposition[9][9];
unsigned propositionCandidate0[9][9];
unsigned propositionCandidate1[9][9];

bool vectorsValidate = false;

GLint modelingMatrixLoc[3];
GLint viewingMatrixLoc[3];
GLint projectionMatrixLoc[3];
GLint eyePosLoc[3];
GLint lightPosLoc[3];
GLint kdLoc[3];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix = glm::translate(glm::mat4(1.f), glm::vec3(-0.5, -0.5, -0.5));

glm::vec3 eyePos = glm::vec3(0, 20, 24);
glm::vec3 lightPos = glm::vec3(0, 10, 24);

glm::vec3 kdGround(0.334, 0.288, 0.635); // this is the ground color in the demo
glm::vec3 kdCubes(0.86, 0.11, 0.31);

int activeProgramIndex = 0;

// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

// For reading GLSL files
bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

GLuint createVS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

	return fs;
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[3]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 0, &face))
    //if (FT_New_Face(ft, "/usr/share/fonts/truetype/gentium-basic/GenBkBasR.ttf", 0, &face)) // you can use different fonts
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            (GLuint) face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();
	gProgram[2] = glCreateProgram();
    gProgram[3] = glCreateProgram();
    
	// Create the shaders for both programs

    GLuint vs0 = createVS("vert_ground.glsl"); // for cube at ground shading
    GLuint fs0 = createFS("frag_ground.glsl");

	GLuint vs1 = createVS("vert2.glsl"); // for border shading
	GLuint fs1 = createFS("frag2.glsl");

	GLuint vs2 = createVS("vert_element.glsl");  // for elements
	GLuint fs2 = createFS("frag_element.glsl");

    GLuint vs3 = createVS("vert_text.glsl");  // for text shading
	GLuint fs3 = createFS("frag_text.glsl");

	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs0);
	glAttachShader(gProgram[0], fs0);

	glAttachShader(gProgram[1], vs1);
	glAttachShader(gProgram[1], fs1);

	glAttachShader(gProgram[2], vs2);
	glAttachShader(gProgram[2], fs2);

    glAttachShader(gProgram[3], vs3);
	glAttachShader(gProgram[3], fs3);

	// Link the programs

    for (int i = 0; i < 4; ++i)
    {
        glLinkProgram(gProgram[i]);
        GLint status;
        glGetProgramiv(gProgram[i], GL_LINK_STATUS, &status);

        if (status != GL_TRUE)
        {
            cout << "Program link failed: " << i << endl;
            exit(-1);
        }
    }


	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 3; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
		lightPosLoc[i] = glGetUniformLocation(gProgram[i], "lightPos");
		kdLoc[i] = glGetUniformLocation(gProgram[i], "kd");

        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(modelingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
        glUniform3fv(eyePosLoc[i], 1, glm::value_ptr(eyePos));
        glUniform3fv(lightPosLoc[i], 1, glm::value_ptr(lightPos));
        glUniform3fv(kdLoc[i], 1, glm::value_ptr(kdCubes));
	}
}


std::vector<GLuint> triangleIndices, lineIndices;
bool initialized = false;
GLuint vao, vertexAttribBuffer, indexBuffer;
std::vector<GLfloat> vertexData, normalData;

int pi_triangleIndices_size;
int pi_lineIndices_size;
int pi_vertexData_size;
int pi_normalData_size;




bool getNthBit(unsigned int num, unsigned int n) {
    return (num & (1u << n)) >> n;
}




void addBox(
    glm::vec3 blf, glm::vec3 brf, glm::vec3 trf, glm::vec3 tlf,
    glm::vec3 blb, glm::vec3 brb, glm::vec3 trb, glm::vec3 tlb
);
void deleteNonGroundElements();

void validateVectors(){
    if(!vectorsValidate){
        deleteNonGroundElements();
        for(int i = 0; i< 9 ; i++){
            for(int j = 0; j< 9 ; j++){
                unsigned stateSignifier = state[i][j];
                unsigned propositionSignifier = proposition[i][j];

                for(int l = 0; l < 12; l++){
                    unsigned presence = getNthBit(stateSignifier,l) || getNthBit(propositionSignifier,l);
                    if(presence){
                        glm::vec3 center(-4.0f+(float)i,0.5f+(float)l,-4.0f+(float)j);
                        addBox(
                            glm::vec3(center[0]-0.5f,center[1]-0.5f,center[2]+0.5f),
                            glm::vec3(center[0]+0.5f,center[1]-0.5f,center[2]+0.5f),
                            glm::vec3(center[0]+0.5f,center[1]+0.5f,center[2]+0.5f),
                            glm::vec3(center[0]-0.5f,center[1]+0.5f,center[2]+0.5f),
                            glm::vec3(center[0]-0.5f,center[1]-0.5f,center[2]-0.5f),
                            glm::vec3(center[0]+0.5f,center[1]-0.5f,center[2]-0.5f),
                            glm::vec3(center[0]+0.5f,center[1]+0.5f,center[2]-0.5f),
                            glm::vec3(center[0]-0.5f,center[1]+0.5f,center[2]-0.5f)
                        );
                    }
                }
            }
        }
        vectorsValidate = true;
    }
}



bool buffersValidate = false;

void validateBuffers(){
    validateVectors();
    

    if(! buffersValidate){    
           

        // Update buffers
        glBindBuffer(GL_ARRAY_BUFFER, vertexAttribBuffer);//sus
        glBufferData(GL_ARRAY_BUFFER, 
                    (vertexData.size() + normalData.size()) * sizeof(GLfloat), 
                    nullptr, GL_STATIC_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, 
                        vertexData.size() * sizeof(GLfloat), vertexData.data());
        glBufferSubData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(GLfloat),
                        normalData.size() * sizeof(GLfloat), normalData.data());

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
                    (triangleIndices.size() + lineIndices.size()) * sizeof(GLuint), 
                    nullptr, GL_STATIC_DRAW);
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
                        triangleIndices.size() * sizeof(GLuint), triangleIndices.data());
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, triangleIndices.size() * sizeof(GLuint),
                        lineIndices.size() * sizeof(GLuint), lineIndices.data());

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 
                            (void*)(vertexData.size() * sizeof(GLfloat)));
        buffersValidate = true;
    }
}
void addBox(
    glm::vec3 blf, glm::vec3 brf, glm::vec3 trf, glm::vec3 tlf,
    glm::vec3 blb, glm::vec3 brb, glm::vec3 trb, glm::vec3 tlb
) {
    if (!initialized) {
                glGenVertexArrays(1, &vao);
                glBindVertexArray(vao);

                glEnableVertexAttribArray(0);
                glEnableVertexAttribArray(1);

                glGenBuffers(1, &vertexAttribBuffer);
                glGenBuffers(1, &indexBuffer);

                glBindBuffer(GL_ARRAY_BUFFER, vertexAttribBuffer);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

                initialized = true;
            }
    // Append 24 vertices (4 for each face)
    std::vector<glm::vec3> newVertices = {
        // Front face
        blf, brf, trf, tlf,
        // Back face
        brb, blb, tlb, trb,
        // Left face
        blb, blf, tlf, tlb,
        // Right face
        brf, brb, trb, trf,
        // Top face
        tlf, trf, trb, tlb,
        // Bottom face
        blf, blb, brb, brf
    };


    for (int i=0; i<newVertices.size(); i++) {
        vertexData.push_back(newVertices[i].x);
        vertexData.push_back(newVertices[i].y);
        vertexData.push_back(newVertices[i].z);
    }

    // Append normal data for each vertex
    std::vector<glm::vec3> normals = {
        // Front face
        glm::vec3(0.0f, 0.0f, 1.0f),
        glm::vec3(0.0f, 0.0f, 1.0f),
        glm::vec3(0.0f, 0.0f, 1.0f),
        glm::vec3(0.0f, 0.0f, 1.0f),
        // Back face
        glm::vec3(0.0f, 0.0f, -1.0f),
        glm::vec3(0.0f, 0.0f, -1.0f),
        glm::vec3(0.0f, 0.0f, -1.0f),
        glm::vec3(0.0f, 0.0f, -1.0f),
        // Left face
        glm::vec3(-1.0f, 0.0f, 0.0f),
        glm::vec3(-1.0f, 0.0f, 0.0f),
        glm::vec3(-1.0f, 0.0f, 0.0f),
        glm::vec3(-1.0f, 0.0f, 0.0f),
        // Right face
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(1.0f, 0.0f, 0.0f),
        // Top face
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        // Bottom face
        glm::vec3(0.0f, -1.0f, 0.0f),
        glm::vec3(0.0f, -1.0f, 0.0f),
        glm::vec3(0.0f, -1.0f, 0.0f),
        glm::vec3(0.0f, -1.0f, 0.0f)
    };

    for (int i=0; i<normals.size(); i++) {
        normalData.push_back(normals[i].x);
        normalData.push_back(normals[i].y);
        normalData.push_back(normals[i].z);
    }

    // Calculate new indices
    GLuint baseIndex = vertexData.size()/3-24;
    std::vector<GLuint> newTriangleIndices;
    for (int i = 0; i < 6; ++i) { // 6 faces
        GLuint startIdx = baseIndex + i * 4 ;
        triangleIndices.push_back(startIdx);
        triangleIndices.push_back(startIdx + 1);
        triangleIndices.push_back(startIdx + 2);

        triangleIndices.push_back(startIdx);
        triangleIndices.push_back(startIdx + 2);
        triangleIndices.push_back(startIdx + 3);
    }


    // Line indices
    std::vector<GLuint> newLineIndices = {
        7, 6, 6, 5,
        4, 5, 4, 7,
        2, 3, 3, 0,
        0, 1, 1, 2,
        0, 5, 3, 6,
        2, 7, 1, 4
    };
    for (GLuint idx : newLineIndices) {
        lineIndices.push_back(baseIndex + idx);
    }
    buffersValidate = false;
}



void initGround(){
    for(int i = -4 ; i < 5 ; i++){    
        for(int j = -4 ; j < 5 ; j++){
            glm::vec3 center(0.0f+(float)i,-0.25f,0.0f+(float)j);
            addBox(
                glm::vec3(center[0]-0.5f,center[1]-0.25f,center[2]+0.5f),
                glm::vec3(center[0]+0.5f,center[1]-0.25f,center[2]+0.5f),
                glm::vec3(center[0]+0.5f,center[1]+0.25f,center[2]+0.5f),
                glm::vec3(center[0]-0.5f,center[1]+0.25f,center[2]+0.5f),
                glm::vec3(center[0]-0.5f,center[1]-0.25f,center[2]-0.5f),
                glm::vec3(center[0]+0.5f,center[1]-0.25f,center[2]-0.5f),
                glm::vec3(center[0]+0.5f,center[1]+0.25f,center[2]-0.5f),
                glm::vec3(center[0]-0.5f,center[1]+0.25f,center[2]-0.5f)
            );
        }
    }
    
    pi_triangleIndices_size = triangleIndices.size();
    pi_lineIndices_size     = lineIndices.size();
    pi_vertexData_size      = vertexData.size();
    pi_normalData_size      = normalData.size();

}

void deleteNonGroundElements(){
    while(pi_triangleIndices_size != triangleIndices.size()){
        triangleIndices.pop_back();
    }

    while(pi_lineIndices_size != lineIndices.size()){
        lineIndices.pop_back();
    }

    while(pi_vertexData_size != vertexData.size()){
        vertexData.pop_back();
    }

    while(pi_normalData_size != normalData.size()){
        normalData.pop_back();
    }

    buffersValidate = false;

}

void assignment9by9(unsigned assignee[9][9], unsigned assigner[9][9]){
    for(int i = 0 ; i < 9 ; i++){
        for(int j = 0 ; j < 9 ; j++){
            assignee[i][j] = assigner[i][j];
        }
    }
}

bool gameEnded = false;
void endGame(){
    gameEnded = true;
}

void randomSetupProposition(){
      
    srand(time(NULL));
    if(rand()%2){
        assignment9by9(proposition,propositionCandidate0);
    }else{
        assignment9by9(proposition,propositionCandidate1);
    }
    for(int i = 0 ; i <9; i++){
        for(int j = 0 ; j <9; j++){
            if(proposition[i][j] & state[i][j]){
                endGame();
            }
        }
    }
    
    vectorsValidate = false;
}
    

void init() 
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    // polygon offset is used to prevent z-fighting between the cube and its borders
    glPolygonOffset(0.5, 0.5);
    glEnable(GL_POLYGON_OFFSET_FILL);

    initShaders();

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            state[i][j] = 0;
            proposition[i][j] = 0;
            propositionCandidate0[i][j] = 0;
            propositionCandidate1[i][j] = 0;
            if( i >= 3 && i <= 5){
                if( j >= 3 && j <= 5){
                    propositionCandidate0[i][j] = 7<<9;
                    propositionCandidate1[i][j] = 0b001111 << 6;
                }
            }
        }
    }
    randomSetupProposition();

    initGround();
    initFonts(gWidth, gHeight);
}

void drawBoxesGround()
{
	glUseProgram(gProgram[0]);
    glDrawElements(GL_TRIANGLES, pi_triangleIndices_size, GL_UNSIGNED_INT, 0);
}

void drawCubes(){
    glUseProgram(gProgram[2]);
    glDrawElements(GL_TRIANGLES, triangleIndices.size()-pi_triangleIndices_size, GL_UNSIGNED_INT, (void*)(pi_triangleIndices_size*sizeof(GLuint)));
}

void drawBoxesEdges()
{
    glLineWidth(3);

	glUseProgram(gProgram[1]);
    glDrawElements(GL_LINES, lineIndices.size(), GL_UNSIGNED_INT, 
               (void*)(triangleIndices.size() * sizeof(GLuint)));
    
}

void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state	
    glUseProgram(gProgram[3]);
    glUniform3f(glGetUniformLocation(gProgram[3], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}


int points = 0;
void display()
{
    validateBuffers();
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    drawBoxesGround();
    drawCubes();
    drawBoxesEdges();
    
    string topLeftText = sideNames[sideIndex];
    string topRightText = "Points:"+to_string(points);
    string gameEndText = "Game Over";

    float scale = 0.75;
    renderText(topLeftText, 0, 1000-30 , scale, glm::vec3(1, 1, 0));
    renderText(topRightText, 600-130, 1000-30 , scale, glm::vec3(1, 1, 0));

    if(gameEnded){
        renderText(gameEndText, 600/2-300, 1000/2-30 , scale*3, glm::vec3(1, 1, 0));        
    }

    assert(glGetError() == GL_NO_ERROR);
}

void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

	// Use perspective projection

	float fovyRad = (float) (45.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, gWidth / (float) gHeight, 1.0f, 100.0f);

    // always look toward (0, 0, 0)
	viewingMatrix = glm::lookAt(eyePos, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));

    for (int i = 0; i < 3; ++i)
    {
        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(projectionMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(viewingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
    }
}

glm::vec3 eyePositionCompute(float rotationAngle){
    float rotationRadians = rotationAngle/180*M_PI;
    glm::mat4 M = glm::rotate(glm::mat4(1.0f), rotationRadians, glm::vec3(0, 1, 0));
    glm::vec4 temp = M * glm::vec4(eyePos,1.0f);
    return glm::vec3(temp);
}

glm::vec3 relativeLs(0,-5,0);

void handleEyePositionChange(GLFWwindow* window){
    viewingMatrix = glm::lookAt(eyePos, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
    for (int i = 0; i < 3; ++i)
    {
        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(projectionMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(viewingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
        glUniform3f(eyePosLoc[i],eyePos[0],eyePos[1],eyePos[2]);
        glUniform3f(lightPosLoc[i],eyePos[0]+relativeLs[0],eyePos[1]+relativeLs[1],eyePos[2]+relativeLs[2]);
    }
    display();
    glfwSwapBuffers(window);
}

unsigned postMovementProposition[9][9];
void tryUnitMovementHorizontal(bool positive = false){
    for(int j = 0 ; j < 9 ; j++){
        if(positive){
            if(proposition[8][j]){
                return;
            }
        }else{    
            if(proposition[0][j]){
                return;
            }
        }
    }
    assignment9by9(postMovementProposition,proposition);
    if(positive){
        for(int i = 8 ; i >= 1 ; i--){
            for(int j = 0; j < 9 ; j++){
                postMovementProposition[i][j] = postMovementProposition[i-1][j]; 
            }
        }
        for(int j = 0; j<9; j++){
            postMovementProposition[0][j] = 0;
        }
    }else{
        for(int i = 0 ; i < 8 ; i++){
            for(int j = 0; j < 9 ; j++){
                postMovementProposition[i][j] = postMovementProposition[i+1][j]; 
            }
        }
        for(int j = 0; j<9; j++){
            postMovementProposition[8][j] = 0;
        }
    }

    bool noColusion = true;
    for(int i = 0 ; i < 9 ; i++){
        for(int j = 0 ; j < 9 ; j++){
            if(state[i][j] & postMovementProposition[i][j]){
                noColusion = false;
                break;
            }
        }
        if(! noColusion){
            break;
        }
    }
    
    if(noColusion){
        assignment9by9(proposition,postMovementProposition);
        vectorsValidate = false;
    }

}

void tryUnitMovementVertical(bool positive = false){
    for(int i = 0 ; i < 9 ; i++){
        if(positive){
            if(proposition[i][8]){
                return;
            }
        }else{    
            if(proposition[i][0]){
                return;
            }
        }
    }
    assignment9by9(postMovementProposition,proposition);
    if(positive){
        for(int j = 8 ; j >= 1 ; j--){
            for(int i = 0; i < 9 ; i++){
                postMovementProposition[i][j] = postMovementProposition[i][j-1]; 
            }
        }
        for(int i = 0; i<9; i++){
            postMovementProposition[i][0] = 0;
        }
    }else{
        for(int j = 0 ; j < 8 ; j++){
            for(int i = 0; i < 9 ; i++){
                postMovementProposition[i][j] = postMovementProposition[i][j+1]; 
            }
        }
        for(int i = 0; i<9; i++){
            postMovementProposition[i][8] = 0;
        }
    }

    bool noColusion = true;
    for(int j = 0 ; j < 9 ; j++){
        for(int i = 0 ; i < 9 ; i++){
            if(state[i][j]&postMovementProposition[i][j]){
                noColusion = false;
                break;
            }
        }
        if(! noColusion){
            break;
        }
    }
    
    if(noColusion){
        assignment9by9(proposition,postMovementProposition);
        vectorsValidate = false;
    }
}


void handlePostFallStateClearence(){
    

}

unsigned postFallProposition[9][9];
void initiateFall(){
    assignment9by9(postFallProposition,proposition);
    bool fallLegal = true;
    for(int i = 0 ; i <9 ; i++){
        for(int j = 0 ; j <9 ; j++){
            if(postFallProposition[i][j]&1){
                fallLegal = false;
                break;
            }
        }
        if(!fallLegal){
            break;
        }
    }

    for(int i = 0 ; i <9 ; i++){
        for(int j = 0 ; j <9 ; j++){
            postFallProposition[i][j] = postFallProposition[i][j] >> 1; 
        }
    }

    if(fallLegal){
        for(int i = 0 ; i <9 ; i++){
            for(int j = 0 ; j <9 ; j++){
                if(postFallProposition[i][j] & state[i][j]){
                    
                    fallLegal = false;
                    break;
                }
            }
            if(!fallLegal){
                break;
            }
        }
    }


    if(!fallLegal){
        
        for(int i = 0 ; i <9 ; i++){
            for(int j = 0 ; j <9 ; j++){
                state [i][j] = state[i][j] | proposition[i][j]; 
                vectorsValidate = false;
            }
        }
        //handlePostFallStateClearence();
        randomSetupProposition();
    }else{
        assignment9by9(proposition, postFallProposition);
        vectorsValidate = false;
    }


}
int counter_maximum = 120;
int counter = 120;
int speed_maximum = 5;
int speed   = 1;


void changeSpeed(bool positive){
    if(positive){
        if(speed < speed_maximum){
            speed += 1;
        }
    }else{
        if(speed > 0){
            speed -= 1;
        }
    }
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_W && action == GLFW_PRESS)
    {
        if(gameEnded){
            return;
        }
        changeSpeed(false);
    }

    if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        if(gameEnded){
            return;
        }
        changeSpeed(true);
    }

    if ((key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE) && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    
    if (key == GLFW_KEY_H && action == GLFW_PRESS)
    {   
        glm::vec3 eventual = eyePositionCompute(-90);

        for(int i = 0; i < turningGraduality; i++){
            eyePos = eyePositionCompute(-90.0f/turningGraduality);
            handleEyePositionChange(window);
        }
        eyePos = eventual;
        handleEyePositionChange(window);
        sideIndex = (sideIndex+3)%4;  
    }
    
    if (key == GLFW_KEY_K && action == GLFW_PRESS)
    {
        glm::vec3 eventual = eyePositionCompute(90);

        for(int i = 0; i < turningGraduality; i++){
            eyePos = eyePositionCompute(90.0f/turningGraduality);
            handleEyePositionChange(window);
        }
        eyePos = eventual;
        handleEyePositionChange(window);   
        sideIndex = (sideIndex+1)%4;  
    }

    if (key == GLFW_KEY_A && action == GLFW_PRESS)
    {
        if(gameEnded){
            return;
        }
        switch (sideIndex) {
            case 0:
                tryUnitMovementHorizontal(false);
                break;
            case 1:
                tryUnitMovementVertical(true);
                break;
            case 2:
                tryUnitMovementHorizontal(true);
                break;
            case 3:
                tryUnitMovementVertical(false);
                break;
            default:
                break;
            }
        }
    
    if (key == GLFW_KEY_D && action == GLFW_PRESS)
    {
        if(gameEnded){
            return;
        }
        switch (sideIndex) {
            case 0:
                tryUnitMovementHorizontal(true);
                break;
            case 1:
                tryUnitMovementVertical(false);
                break;
            case 2:
                tryUnitMovementHorizontal(false);
                break;
            case 3:
                tryUnitMovementVertical(true);
                break;
            default:
                break;
        }
        
    }


}

void mainLoop(GLFWwindow* window)
{
    while (!glfwWindowShouldClose(window))
    {   
        if(! gameEnded){
            counter -= speed;
            if(counter < 0){
                counter = counter_maximum;
                initiateFall();
            }
        }
        

        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(gWidth, gHeight, "tetrisGL", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, gWidth, gHeight); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
