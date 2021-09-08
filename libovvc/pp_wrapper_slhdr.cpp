/**
 * \example SLHDRDecode.cpp
 */

#include <cstdio>

#include "SLHDRPostprocessor/SLPostprocessorContext.h"
#include "SLHDRCommon/SLCommon.h"

extern "C" void 
pp_init_slhdr_lib(void** pslhdr_context)
{
    //Request PQ as transfer function for reconstructed HDR
    SLHDR::SLPostprocessorContextProperties postprocessorProperties;
    postprocessorProperties.HDRTransferFunction = SLHDR::PQ;

    //print postprocessor properties
    postprocessorProperties.print();

    //create a the cpu list to specify on which cpu threads from the thread pool will be binded (we use 4 cpu here)
    std::vector<u8> cpulist;
    for (u8 cpu = 0; cpu < SLHDR::getMaxThreadAvailableCount(); cpu++) 
        cpulist.push_back(cpu);

    //instanciate postprocessor context with filled properties
    SLHDR::SLPostprocessorContext* dsl = new SLHDR::SLPostprocessorContext(postprocessorProperties,cpulist);
    *pslhdr_context = dsl;

    //we set the postprocessor display adaptation profile to the medium level, low and high are also available
    dsl->setPostprocessorProfile(SLHDR::SLPostprocessorProfile::med);

}

extern "C" void 
pp_uninit_slhdr_lib(void* slhdr_context)
{
    SLHDR::SLPostprocessorContext* dsl = (SLHDR::SLPostprocessorContext*)slhdr_context;
    delete dsl;
}

extern "C" void 
pp_slhdr_no_filter(void* slhdr_context, int16_t** sdr_pic, int16_t** hdr_pic, 
                uint8_t* SEIPayload, int pic_width, int pic_height)
{
    return;
}

extern "C" void 
pp_sdr_to_hdr(void* slhdr_context, int16_t** sdr_pic, int16_t** hdr_pic, 
                uint8_t* SEIPayload, int pic_width, int pic_height)
{

    SLHDR::SLPostprocessorContext* dsl = (SLHDR::SLPostprocessorContext*)slhdr_context;

    //allocate memory for the HD YUV SDR 420 input
    // SLHDR::SLBuffer inPic(1920, 1080,2,SLHDR::p420);
    SLHDR::SLBuffer inPic(sdr_pic[0], sdr_pic[1], sdr_pic[2], pic_width, pic_height, 2, SLHDR::p420);

    //allocate memory for the reconstructed HD YUV HDR v210 output
    // SLHDR::SLBuffer reconstructedHDRPic(1920, 1080,2, SLHDR::p420);
    SLHDR::SLBuffer reconstructedHDRPic(hdr_pic[0], hdr_pic[1], hdr_pic[2], pic_width, pic_height, 2, SLHDR::p420);

    SLHDR::SLError err = dsl->decode(inPic, SEIPayload, reconstructedHDRPic);

    if(err)
    {
        printf("postprocessor error : %s\n",dsl->getErrorString().c_str());
        return ;
    }

    //print input stream characteristics. only works after one frame has been decoded.
    // dsl.getStreamProperties().print();

 
    printf("done !\n");

    return ;
}