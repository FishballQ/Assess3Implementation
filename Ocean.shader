Shader "Ocean"
{
    Properties
    {
        _OceanColorShallow ("Ocean Color Shallow", Color) = (1, 1, 1, 1) // shallow color
        _OceanColorDeep ("Ocean Color Deep", Color) = (1, 1, 1, 1)       // deep color
        _BubblesColor ("Bubbles Color", Color) = (1, 1, 1, 1)            // bubbles color
        _Specular ("Specular", Color) = (1, 1, 1, 1)                     // highlight
        _Gloss ("Gloss", Range(8.0, 256)) = 20                           // gloss
        _FresnelScale ("Fresnel Scale", Range(0, 1)) = 0.5               // fresnel effect
        _Roughness ("roughness", Range(0, 50)) = 0                       // reflection roughness
        _Displace ("Displace", 2D) = "black" { }                         // offset texture
        _Normal ("Normal", 2D) = "black" { }                             // normal texture
        _Bubbles ("Bubbles", 2D) = "black" { }                           // bubbles texture
        _DepthMaxDistance ("Depth Max Distance", float) = 150.0          // depth
        _DepthVisibility ("Depth Visibility", float) = 100.0             // depth visibility
        _FoamDistance("Foam Distance", Float) = 100                      // foam detection
        //replace last foam distance variable
        _FoamMaxDistance("Foam Maximum Distance", Float) = 0.4
        _FoamMinDistance("Foam Minimum Distance", Float) = 0.04
        _OceanColorSSS ("Ocean Color SSS", Color) = (1, 1, 1, 1)
        _SSSPow ("SSS Pow", Range(.001, 5)) = .5
        _SSSHeight ("SSS Height", float) = .5
        _SSSIntensity ("SSS Intensity", float) = .5
        _EmissionStrength ("Emission Strength", Range(0.001, 100)) = .5
    }
    SubShader
    {
        Tags { "RenderType" = "Transparent" "LightMode" = "ForwardBase" "Queue" = "Transparent" }
        LOD 100
        
        Pass
        {
            Blend SrcAlpha OneMinusSrcAlpha 
            ZWrite Off 
            
            CGPROGRAM
            
            #pragma vertex vert
            #pragma fragment frag
            #pragma enable_d3d11_debug_symbols
            
            #include "UnityCG.cginc"
            #include "Lighting.cginc"
            
            struct appdata
            {
                float4 vertex: POSITION;
                float2 uv: TEXCOORD0;
            };
            
            struct v2f
            {
                float4 pos: SV_POSITION;
                float2 uv: TEXCOORD0;
                float3 worldPos: TEXCOORD1;
                float4 screenPosition : TEXCOORD2;
            };
            
            fixed4 _OceanColorShallow;
            fixed4 _OceanColorDeep;
            fixed4 _BubblesColor;
            fixed4 _Specular;
            float _Gloss;
            float _Roughness;
            fixed _FresnelScale;
            sampler2D _Displace;
            sampler2D _Normal;
            sampler2D _Bubbles;
            float4 _Displace_ST;
            sampler2D _CameraDepthTexture;
            fixed _DepthMaxDistance;
            float _FoamDistance;
            sampler2D _CameraNormalsTexture;
            float _FoamMaxDistance;
            float _FoamMinDistance;
            sampler2D _CameraOpaqueTexture;
            sampler2D sampler_CameraOpaqueTexture;
            sampler2D _WaveResult;
            sampler2D _Tex;
            float _DepthVisibility;
            float _SSSPow;
            float _SSSIntensity;
            float _SSSHeight;
            float _EmissionStrength;
            fixed4 _OceanColorSSS;
            
            v2f vert(appdata v)
            {
                v2f o;
                // uv compares with the tiling and offset value of material, to ensure the correct tiling and offset
                o.uv = TRANSFORM_TEX(v.uv, _Displace);
                // x, y are texture coordinate in uv space, w is mip sample, 0 is the highest resolution
                float4 displace = tex2Dlod(_Displace, float4(o.uv, 0, 0));
                
                // get mesh coordinates
                float4 waveTransmit = tex2Dlod(_WaveResult, float4(v.uv, 0, 0));
                float waveHeight = DecodeFloatRGBA(waveTransmit);
                float4 x = tex2Dlod(_Tex, float4(v.uv, 0, 0));

                v.vertex.xyz += float4(displace.xyz, 0);
                v.vertex.y += waveHeight * 5;
                if(x.r > 0)
                    v.vertex.y = 0;

                o.pos = UnityObjectToClipPos(v.vertex);
                // world
                o.worldPos = mul(unity_ObjectToWorld, v.vertex).xyz;
                // screen 
                o.screenPosition = ComputeScreenPos(o.pos);
                return o;
            }

            fixed meanFresnel(float cosThetaV, float sigmaV)
            {
                return pow(1.0 - cosThetaV, 5.0 * exp(-2.69 * sigmaV)) / (1.0 + 22.7 * pow(sigmaV, 1.5));
            }

            // V, N in world space
            float meanFresnel(fixed3 V, fixed3 N, fixed2 sigmaSq) {
                fixed length2 = 1.0 - V.z * V.z;
                fixed cosPhi2 = V.x * V.x / length2;
                fixed sinPhi2 = V.y * V.y / length2;
                fixed2 t = fixed2(cosPhi2, sinPhi2); // cos^2 and sin^2 of view direction
                float sigmaV2 = dot(t, sigmaSq);     // slope variance in view direction
                return meanFresnel(dot(V, N), sqrt(sigmaV2));
            }

            // V, N, Tx, Ty in world space
            fixed3 U(fixed2 zeta, fixed3 V, fixed3 N) {
                fixed3 f = normalize(fixed3(-zeta, 1.0)); // tangent space
                fixed3 F = f.x + f.y + f.z * N; // world space
                fixed3 R = 2.0 * dot(F, V) * F - V;
                return R;
            }

            // V, N, Tx, Ty in world space;
            fixed3 meanSkyRadiance(fixed3 V, fixed3 N, fixed2 sigmaSq) {
                fixed3 u0 = U(fixed2(0.0, 0.0), V, N);
                float roughness = _Roughness * (1.7 - 0.7 * _Roughness);           //roughness
                half mip = roughness * 6;
                half4 rgbm = UNITY_SAMPLE_TEXCUBE_LOD(unity_SpecCube0, u0, mip);   //reflection probes
                half3 sky = DecodeHDR(rgbm, unity_SpecCube0_HDR);                  //ambient light reflection
                return sky;
            }

            // assumes x>0
            float erfc(float x) {
                return 2.0 * exp(-x * x) / (2.319 * x + sqrt(4.0 + 1.52 * x * x));
            }

            float Lambda(float cosTheta, float sigmaSq) {
                float v = cosTheta / sqrt((1.0 - cosTheta * cosTheta) * (2.0 * sigmaSq));
                return max(0.0, (exp(-v * v) - v * sqrt(UNITY_PI) * erfc(v)) / (2.0 * v * sqrt(UNITY_PI)));
                //return (exp(-v * v)) / (2.0 * v * sqrt(M_PI)); // approximate, faster formula
            }

            // L, V, N, Tx, Ty in world space
            float reflectedSunRadiance(fixed3 L, fixed3 V, fixed3 N, fixed2 sigmaSq) {
                fixed3 H = normalize(L + V);

                float zL = dot(L, N); // cos of source zenith angle
                float zV = dot(V, N); // cos of receiver zenith angle
                float zH = dot(H, N); // cos of facet normal zenith angle
                float zH2 = zH * zH;

                sigmaSq = max(sigmaSq, 2e-5);
                float p = exp(-1.0f) / (2.0 * UNITY_PI * sqrt(sigmaSq.x * sigmaSq.y));

                fixed lengthV2 = 1 - V.z * V.z;
                fixed SinV2 = V.y * V.y / lengthV2;
                fixed CosV2 = V.x * V.x / lengthV2;
                float sigmaV2 = sigmaSq.x * CosV2 + sigmaSq.y * SinV2;


                fixed lengthL2 = 1 - L.z * L.z;
                fixed SinL2 = L.y * L.y / lengthL2;
                fixed CosL2 = L.x * L.x / lengthL2;
                float sigmaL2 = sigmaSq.x * CosL2 + sigmaSq.y * SinL2;

                float fresnel = 0.02 + 0.98 * pow(1.0 - dot(V, H), 5.0);

                zL = max(zL, 0.01);
                zV = max(zV, 0.01);

                return fresnel * p / ((1.0 + Lambda(zL, sigmaL2) + Lambda(zV, sigmaV2)) * zV * zH2 * zH2 * 4.0);
                //return p;
            }

            float4 CalculateSSSColor(float3 lightDirection, float3 worldNormal, float3 viewDir, float shadowFactor)
            {
                    float lightStrength = sqrt(saturate(lightDirection.y));
                    float SSSFactor = pow(saturate(dot(viewDir ,lightDirection) )+saturate(dot(worldNormal ,-lightDirection)) ,_SSSPow) * shadowFactor * lightStrength * _EmissionStrength;
                    return _OceanColorSSS * SSSFactor;
            }
            
            fixed4 frag(v2f i): SV_Target
            {
                // depth texture
                float existingDepth01 = tex2Dproj(_CameraDepthTexture, UNITY_PROJ_COORD(i.screenPosition)).r;
                float existingDepthLinear = LinearEyeDepth(existingDepth01);
                // water depth
                float waterDepth = i.screenPosition.w;
                float depthDifference = existingDepthLinear - waterDepth;
                float waterDepthDifference01 = saturate(depthDifference / _DepthMaxDistance);
                // foam depth
	            float foamDepthDifference01 = saturate(depthDifference / _FoamDistance);
                float alpha = saturate(depthDifference / _DepthVisibility);
                alpha = lerp(alpha, 1, foamDepthDifference01);
                
                // get normal from the normal texture, transfer to the world coordinate
                fixed3 normal = UnityObjectToWorldNormal(tex2D(_Normal, i.uv).rgb);

                // UV screen coordinate
                float2 screenUV = i.screenPosition.xy / i.screenPosition.w;

                // calculate the view space normal
                fixed3 waterNormal = UnityObjectToViewPos(normal);
                float3 existingNormal = tex2Dproj(_CameraNormalsTexture, UNITY_PROJ_COORD(i.screenPosition));
                float3 normalDot = saturate(dot(existingNormal, waterNormal));
                // smooth and gradient effect for the foam using Lerp
                float foamDistance = lerp(_FoamMaxDistance, _FoamMinDistance, normalDot);
                
                fixed2 sigmaSq = fixed2(normal.x * normal.x, normal.z * normal.z);
                fixed bubbles = tex2D(_Bubbles, i.uv).r;    // bubble intensity using bubble texture

                fixed3 lightDir = normalize(UnityWorldSpaceLightDir(i.worldPos));  // light direction
                fixed3 viewDir = normalize(UnityWorldSpaceViewDir(i.worldPos));    // view direction
                fixed3 reflectDir = reflect(-viewDir, normal);                     //reflect direction
                reflectDir = BoxProjectedCubemapDirection(reflectDir, i.worldPos, unity_SpecCube0_ProbePosition, unity_SpecCube0_BoxMin, unity_SpecCube0_BoxMax);

                float roughness = _Roughness * (1.7 - 0.7 * _Roughness);                  // roughness
                half mip = roughness * 6;
                half4 rgbm = UNITY_SAMPLE_TEXCUBE_LOD(unity_SpecCube0, reflectDir, mip);  //reflection probes
                half3 sky = DecodeHDR(rgbm, unity_SpecCube0_HDR);                         // ambient light reflection
                // fresnel
                fixed fresnel =  0.02 + 0.98 * meanFresnel(viewDir, normal, sigmaSq);
                fresnel = saturate(_FresnelScale + (1 - _FresnelScale) * pow(1 - dot(normal, viewDir), 5));
                // reflective direction, blend shallow color and deep color based on the view direction
                half facing = saturate(dot(viewDir, normal));                
                fixed3 oceanColor = lerp(_OceanColorShallow, _OceanColorDeep, facing);

                fixed3 ambient = UNITY_LIGHTMODEL_AMBIENT.rgb;

                // bubbles color using diffuse 
                fixed3 bubblesDiffuse = _BubblesColor.rbg * _LightColor0.rgb * saturate(dot(lightDir, normal));
                // ocean color using semi vectorial, reflected light and highlights
                fixed3 halfDir = normalize(lightDir + viewDir);
                fixed3 specular = _LightColor0.rgb * _Specular.rgb * pow(max(0, dot(normal, halfDir)), _Gloss);
                // ocean diffuse 
                fixed3 oceanDiffuse = oceanColor * _LightColor0.rgb * saturate(dot(lightDir, normal));
                fixed3 diffuse = lerp(oceanDiffuse, bubblesDiffuse, bubbles);

                // final ocean color: ambient light, diffuse, reflection by sky (according to fresnel)
                fixed3 col = fixed3(0.0f, 0.0f, 0.0f);
                // sun color
                fixed3 Lsun = _LightColor0.rgb / 10000.0f;
 
                col += reflectedSunRadiance(lightDir, viewDir, normal, sigmaSq) * Lsun;
                // Sky color, light
                col += fresnel * meanSkyRadiance(viewDir, normal, sigmaSq);
                // water color, Refracted light
                fixed3 Lsea = diffuse * sky;
                col += Lsea * (1 - fresnel);

                col = lerp(col, bubblesDiffuse.rbg, bubbles);
                // ocean ambient + highlight
                col = ambient + lerp(diffuse, sky, fresnel) + specular;

                fixed4 sss = CalculateSSSColor(lightDir, normal, viewDir,_SSSIntensity);
                col = (1 - sss.a) * col + sss.a * sss;
                
                //col = fixed3(depthDifference, depthDifference, depthDifference);
                float3 waterColor = lerp(col, col * 0.1 + _OceanColorDeep, waterDepthDifference01);
                float3 foamColor = lerp(_BubblesColor.rbg, waterColor, foamDepthDifference01);
                //col = waterColor;
                float4 waveTransmit = tex2Dlod(_WaveResult, float4(i.uv, 0, 0));
                float waveHeight = DecodeFloatRGBA(waveTransmit);

                col = foamColor;

                if(waveHeight < 0.1)
                    return fixed4(col, alpha);
                else
                {
                    sky = lerp(col, _BubblesColor, waveHeight / 2);
                    return fixed4(sky, alpha);
                }
                return fixed4(col, alpha);

                //col = foamColor;

                //return fixed4(col, alpha);
            }
            
            ENDCG
            
        }
    }
}
