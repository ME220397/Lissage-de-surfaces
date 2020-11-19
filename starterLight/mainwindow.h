#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value; Color faceShadingColor;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;


enum DisplayMode {Normal, TemperatureMap, ColorShading};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // les 4 fonctions à compléter
    void showSelections(MyMesh* _mesh);
    void showSelectionsNeighborhood(MyMesh* _mesh);
    void showPath(MyMesh* _mesh, int v1, int v2);
    void showBorder(MyMesh* _mesh);

    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    // Fonctions TPs
    void operateur_laplace_beltrami(MyMesh* _mesh, int choix, double h, double _y);
    MyMesh::Point laplace_beltrami_cot(MyMesh* _mesh, VertexHandle v);
    MyMesh::Point laplace_beltrami_uni(MyMesh* _mesh, VertexHandle v);
    double calcul_poids_cot(MyMesh* _mesh, VertexHandle vi, VertexHandle vj);
    MyMesh::Scalar neighboring_faces_area(MyMesh* _mesh, VertexHandle v);
    MyMesh::Point calc_vector_v_vi(MyMesh* _mesh, VertexHandle v, VertexHandle vi);
    void flou_de_diffusion(MyMesh *_mesh, VertexHandle v, MyMesh::Point vi, double f);
private slots:

    void on_pushButton_chargement_clicked();
    void on_operateur_clicked();

    void on_pushButton_clicked();

    void on_doubleSpinBox_valueChanged(double arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_scale_clicked();

    void on_Unscale_clicked();

    void on_pushButton_2_clicked();

private:

    bool modevoisinage;
    const int COTANGENTE = 0;
    const int UNIFORME = 0;
    MyMesh mesh;
    MyMesh clone;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    double h = -0.001;
    double _y = 0.001;
    bool factor_change = false;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
