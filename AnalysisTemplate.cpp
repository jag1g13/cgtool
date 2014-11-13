#include <iostream>
#include <string>
#include <vector>
#include <gromacs/trajectoryanalysis.h>

using namespace gmx;


//int main(int argc, char **argv) {
//    std::cout << "Hello, world!" << std::endl;
//    return 0;
//}


class AnalysisTemplate : public TrajectoryAnalysisModule
{
    public:
        AnalysisTemplate();
        virtual void initOptions(Options *options,
        TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
        const TopologyInformation &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
        TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();
    private:
        class ModuleData;
        std::string fnDist_;
        double cutoff_;
        Selection refsel_;
        SelectionList sel_;
        AnalysisNeighborhood nb_;
        AnalysisData data_;
        AnalysisDataAverageModulePointer avem_;
};

AnalysisTemplate::AnalysisTemplate()
    : TrajectoryAnalysisModule("template", "Template analysis tool"),
      cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
AnalysisTemplate::initOptions(Options *options,
TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "Gromacs. The advantage of using Gromacs for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by Gromacs. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };
    options->setDescription(desc);
    options->addOption(FileNameOption("o")
                        .filetype(eftPlot).outputFile()
                        .store(&fnDist_).defaultBasename("avedist")
                        .description("Average distances from reference group"));
    options->addOption(SelectionOption("reference")
                        .store(&refsel_).required()
                        .description("Reference group to calculate distances from"));
    options->addOption(SelectionOption("select")
                        .storeVector(&sel_).required().multiValue()
                        .description("Groups to calculate distances to"));
    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                        .description("Cutoff for distance calculation (0 = no cutoff)"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
const TopologyInformation & /*top*/)
{
    nb_.setCutoff(cutoff_);
    data_.setColumnCount(0, sel_.size());
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);
    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        data_.addModule(plotm);
    }
}


void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dh = pdata->dataHandle(data_);
    const Selection &refsel = pdata->parallelSelection(refsel_);
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refsel);
    dh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        const Selection &sel = pdata->parallelSelection(sel_[g]);
        int nr = sel.posCount();
        real frave = 0.0;
        for (int i = 0; i < nr; ++i)
        {
            SelectionPosition p = sel.position(i);
            frave += nbsearch.minimumDistance(p.x());
        }
        frave /= nr;
        dh.setPoint(g, frave);
    }
    dh.finishFrame();
}


void
AnalysisTemplate::finishAnalysis(int /*nframes*/)
{
}


void
AnalysisTemplate::writeOutput()
{
    // We print out the average of the mean distances for each group.
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        fprintf(stderr, "Average mean distance for '%s': %.3f nm\n",
        sel_[g].name(), avem_->average(0, g));
    }
}


/*! \brief
* The main function for the analysis template.
*/

int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(argc, argv);
}