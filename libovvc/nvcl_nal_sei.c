#include <stddef.h>

#include "ovutils.h"
#include "ovmem.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"

enum SEIPayloadtype
{
    BUFFERING_PERIOD                     = 0,
    PICTURE_TIMING                       = 1,
    FILLER_PAYLOAD                       = 3,
    USER_DATA_REGISTERED_ITU_T_T35       = 4,
    USER_DATA_UNREGISTERED               = 5,
    FILM_GRAIN_CHARACTERISTICS           = 19,
    FRAME_PACKING                        = 45,
    PARAMETER_SETS_INCLUSION_INDICATION  = 129,
    DECODING_UNIT_INFO                   = 130,
    DECODED_PICTURE_HASH                 = 132,
    SCALABLE_NESTING                     = 133,
    MASTERING_DISPLAY_COLOUR_VOLUME      = 137,
    DEPENDENT_RAP_INDICATION             = 145,
    EQUIRECTANGULAR_PROJECTION           = 150,
    SPHERE_ROTATION                      = 154,
    REGION_WISE_PACKING                  = 155,
    OMNI_VIEWPORT                        = 156,
    GENERALIZED_CUBEMAP_PROJECTION       = 153,
    FRAME_FIELD_INFO                     = 168,
    SUBPICTURE_LEVEL_INFO                = 203,
    SAMPLE_ASPECT_RATIO_INFO             = 204,
    CONTENT_LIGHT_LEVEL_INFO             = 144,
    ALTERNATIVE_TRANSFER_CHARACTERISTICS = 147,
    AMBIENT_VIEWING_ENVIRONMENT          = 148,
    CONTENT_COLOUR_VOLUME                = 149,
    ANNOTATED_REGIONS                    = 202,
};

struct OVSEIPayload
{
    int type;
    uint32_t size;
};

void
copy_sei_params(OVSEI **dst_p, OVSEI *src)
{   
    if(src){
        if(!(*dst_p))
            *dst_p = ov_mallocz(sizeof(struct OVSEI));

        OVSEI *dst = *dst_p; 
        if(!dst->sei_fg){
            dst->sei_fg = ov_mallocz(sizeof(struct OVSEIFGrain));
        }
        *(dst->sei_fg) =  *(src->sei_fg);
    }
}

/* FIXME find other spec */
struct OVSEIPayload
nvcl_sei_payload(OVNVCLReader *const rdr) {
    int pl_type = 0;
    uint8_t val = 0;
    do
    {
        val = nvcl_read_bits(rdr, 8);
        pl_type += val;
    } while (val==0xFF);

    uint32_t pl_size = 0;
    do
    {
        val = nvcl_read_bits(rdr, 8);
        pl_size += val;
    } while (val==0xFF);

    struct OVSEIPayload payload;
    payload.type = pl_type;
    payload.size = pl_size;
    return payload;
}

void
nvcl_film_grain_read(OVNVCLReader *const rdr, struct OVSEIFGrain *const fg, OVNVCLCtx *const nvcl_ctx)
{
  fg->fg_characteristics_cancel_flag = nvcl_read_flag(rdr);
  if (!fg->fg_characteristics_cancel_flag )
  {
    fg->fg_model_id = nvcl_read_bits(rdr,2);
    fg->fg_separate_colour_description_present_flag = nvcl_read_flag(rdr);
    if (fg->fg_separate_colour_description_present_flag)
    {
        fg->fg_bit_depth_luma_minus8 = nvcl_read_bits(rdr, 3);
        fg->fg_bit_depth_chroma_minus8 = nvcl_read_bits(rdr, 3);
        fg->fg_full_range_flag = nvcl_read_flag(rdr);
        fg->fg_colour_primaries = nvcl_read_bits(rdr, 8);
        fg->fg_transfer_characteristics = nvcl_read_bits(rdr, 8);
        fg->fg_matrix_coeffs = nvcl_read_bits(rdr, 8);
    }
    fg->fg_blending_mode_id = nvcl_read_bits(rdr, 2);
    fg->fg_log2_scale_factor = nvcl_read_bits(rdr, 4);
    for (int c = 0; c<3; c++)
    {
      fg->fg_comp_model_present_flag[c] = nvcl_read_flag(rdr);
    }
    for (int c = 0; c<3; c++)
    {
      if (fg->fg_comp_model_present_flag[c])
      {
        fg->fg_num_intensity_intervals_minus1[c] = nvcl_read_bits(rdr, 8);
        fg->fg_num_model_values_minus1[c] = nvcl_read_bits(rdr, 3);

        for (uint32_t i = 0; i< fg->fg_num_intensity_intervals_minus1[c] + 1 ; i++)
        {
          fg->fg_intensity_interval_lower_bound[c][i] = nvcl_read_bits(rdr, 8);
          fg->fg_intensity_interval_upper_bound[c][i] = nvcl_read_bits(rdr, 8);
          for (uint32_t j = 0; j < fg->fg_num_model_values_minus1[c] + 1; j++)
          {
            fg->fg_comp_model_value[c][i][j] = nvcl_read_s_expgolomb(rdr);
          }
        }
      }
    } // for c
    fg->fg_characteristics_persistence_flag = nvcl_read_flag(rdr);
  } // cancel flag
}

int 
nvcl_decode_nalu_sei(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{   
    if(!nvcl_ctx->sei)
        nvcl_ctx->sei = ov_mallocz(sizeof(struct OVSEI));    
    struct OVSEI* sei = nvcl_ctx->sei;
    
    struct OVSEIPayload payload = nvcl_sei_payload(rdr);
     switch (payload.type)
    {
        uint8_t sei_byte;
        case FILM_GRAIN_CHARACTERISTICS:
            if(!sei->sei_fg)
                sei->sei_fg = ov_mallocz(sizeof(struct OVSEIFGrain));
            nvcl_film_grain_read(rdr, sei->sei_fg, nvcl_ctx);
            break;
        default:
            for (int i = 0; i < payload.size; i++)
            {
                sei_byte = nvcl_read_bits(rdr, 8);
                sei_byte++;
            }
            ov_log(NULL, OVLOG_INFO, "SEI: Unknown prefix message (type = %d) was found!\n", payload.type);
            break;
    }

    return 0;
}
