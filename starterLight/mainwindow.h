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
    void operateur_laplace_beltrami(MyMesh* _mesh, int choix);
    void laplace_beltrami_cot(MyMesh* _mesh, VertexHandle v);
    void laplace_beltrami_uni(MyMesh* _mesh, VertexHandle v);
    MyMesh::Scalar calcul_poids_cot(MyMesh* _mesh, VertexHandle vi, VertexHandle vj);
    MyMesh::Scalar calcul_aire_barycentres(MyMesh* _mesh, VertexHandle v);
    MyMesh::Point direction_v_vi(MyMesh* _mesh, VertexHandle v, VertexHandle vi);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_vertexMoins_clicked();
    void on_pushButton_vertexPlus_clicked();
    void on_pushButton_edgeMoins_clicked();
    void on_pushButton_edgePlus_clicked();
    void on_pushButton_faceMoins_clicked();
    void on_pushButton_facePlus_clicked();
    void on_pushButton_afficherChemin_clicked();
    void on_pushButton_voisinage_clicked();
    void on_pushButton_bordure_clicked();

private:

    bool modevoisinage;
    const int COTANGENTE = 0;
    const int UNIFORME = 0;
    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
