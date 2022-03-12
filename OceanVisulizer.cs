using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using Unity.Collections;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class OceanVisulizer : MonoBehaviour
{
    // water size
    [Header("Size")]
    [Range(3, 11)]
    public int FFTPow = 10;         // power of 2 of the ocean texture size (for FFT) 
    public int MeshSize = 100;		// the number of the ocean mesh size
    public float MeshLength = 512;	// length of the ocean mesh
    private float geometryCellSize; // unit length

    public int FFTSize;			    // the size of the ocean texture, pow(2, FFTPow), Lx = Lz
    
    // mesh components in the ocean object
    private Mesh mesh;              
    private MeshFilter filter;
    private MeshRenderer render;
    
    private Vector3[] vertices;     // vertice positions of the mesh
    private int[] triangles;		// mesh triangles index
    private Vector2[] uvs; 			// uv coordinates

    // water parameters
    [Header("Parameter")]
    public ComputeShader OceanCS;   // shader of computing ocean waves
    
    public float Amplitude = 100;   // amplitude in the phillips spectrum
    [Header("Offset and Bubbles")]
    public float Lambda = 20;              // wave offset
    public float HeightScale = 50;         // wave height
    public float BubblesScale = 1;         // bubble scale
    public float BubblesThreshold = 0.86f; // bubble threshold
    public float TimeScale = 1;            // time
    [Header("Wind")]
    public float WindScale = 2;                    // wind intensity
    public Vector2 Wind = new Vector2(1.0f, 1.0f); // wind speed
    
    // FFT Controllor
    [Range(1, 12)]
    private int ControlM = 6;       // control m , to control FFT transform
    private bool isControlH = false;  // check if calculate horizontal FFT, "false" means use vertical FFT

    // parameters for interactive wave transfer
    [Header("Interactive Wave")]
    private Material m_waterWaveMarkMat;
    private Material m_waveTransmitMat;
    Vector4 m_waveMarkParams;
    private RenderTexture m_waterWaveMarkTexture;
    private RenderTexture m_waveTransmitTexture;
    private RenderTexture m_prevWaveMarkTexture;
    public RenderTexture objectRenderTexture;
    public Camera interactiveCam;
    public Shader addShader;
    private Material addMaterial;
    private RenderTexture TempRT;

    // wave particles
    private Texture2D particlePosTex;
    NativeArray<float> pixData;
    private List<WaveParticle> _waveParticles;
    public Shader waveFilterShader;
    private Material waveFilterMaterial;
    private RenderTexture m_TmpHeightFieldRT;
    private RenderTexture waveParticlePointRT;
    public Material texMaterial;
    
    //material rendering    
    [Header("Material")]
    public Material OceanMaterial;  // ocean rendering material  
    public Material DisplaceXMat;   // offset X material
    public Material DisplaceYMat;   // offset Y material
    public Material DisplaceZMat;   // offset Z material
    public Material DisplaceMat;    // offset
    public Material NormalMat;      // normal material
    public Material BubblesMat;     // bubbles material
    public Material GuassianMat;    // guassian random material
    public Material WaveMat; // interactive wave material

    private int kernelComputeGaussianRandom;            //计算高斯随机数
    private int kernelCreateHeightSpectrum;             //创建高度频谱
    private int kernelCreateDisplaceSpectrum;           //创建偏移频谱
    private int kernelFFTHorizontal;                    //FFT横向
    private int kernelFFTHorizontalEnd;                 //FFT横向，最后阶段
    private int kernelFFTVertical;                      //FFT纵向
    private int kernelFFTVerticalEnd;                   //FFT纵向,最后阶段
    private int kernelTextureGenerationDisplace;        //生成偏移纹理
    private int kernelTextureGenerationNormalBubbles;   //生成法线和泡沫纹理
    private RenderTexture GaussianRandomRT;             //高斯随机数
    private RenderTexture HeightSpectrumRT;             //高度频谱
    private RenderTexture DisplaceXSpectrumRT;          //X偏移频谱
    private RenderTexture DisplaceZSpectrumRT;          //Z偏移频谱
    private RenderTexture DisplaceRT;                   //偏移频谱
    private RenderTexture OutputRT;                     //临时储存输出纹理
    private RenderTexture NormalRT;                     //法线纹理
    private RenderTexture BubblesRT;                    //泡沫纹理
    private float time = 0;

    public Slider a;
    public Slider t;
    public Slider w;
    public Text na;
    public Text nt;
    public Text nw;
    
    private void Awake()
    {
        // set ocean mesh
        GenerateOceanMesh();

        Shader.SetGlobalTexture("_WaveResult", m_waterWaveMarkTexture);
        Shader.SetGlobalFloat("_WaveHeight", WaveHeight);
        interactiveCam.orthographicSize = MeshLength / 2;
        addMaterial = new Material(addShader);
    }
    
    // create ocean mesh using MeshSize and MeshLength 
    private void GenerateOceanMesh()
    {
        // get parameters
        geometryCellSize = MeshLength / MeshSize;
        filter = gameObject.GetComponent<MeshFilter>();
        render = gameObject.GetComponent<MeshRenderer>();
        mesh = new Mesh();
        mesh.name = "Ocean Mesh";
        filter.mesh = mesh;
        render.material = OceanMaterial;
        // initialize mesh 
        vertices = new Vector3[(MeshSize + 1) * (MeshSize + 1)];
        uvs = new Vector2[(MeshSize + 1) * (MeshSize + 1)];
        triangles = new int[MeshSize * MeshSize * 6];
        // set mesh values
        int inx = 0;
        for (int i = 0; i <= MeshSize; i++)
        {
            for (int j = 0; j <= MeshSize; j++)
            {
                int index = i * (MeshSize + 1) + j;
                vertices[index] = new Vector3((j - MeshSize / 2.0f) * geometryCellSize, 0, (i - MeshSize / 2.0f) * geometryCellSize);
                uvs[index] = new Vector2(j / (MeshSize * 1.0f), i / (MeshSize * 1.0f));

                if (i != MeshSize && j != MeshSize)
                {
                    triangles[inx++] = index;
                    triangles[inx++] = index + MeshSize + 1;
                    triangles[inx++] = index + MeshSize + 2;

                    triangles[inx++] = index;
                    triangles[inx++] = index + MeshSize + 2;
                    triangles[inx++] = index + 1;
                }
            }
        }
        // assign for the ocean mesh
        mesh.vertices = vertices;
        mesh.uv = uvs;
        mesh.triangles = triangles;
    }

    private void Start()
    {
        // initialize ComputerShader parameters
        OceanComputerShaderInit();

        m_waterWaveMarkMat = new Material(Shader.Find("Unlit/WaveMarkerShader"));
        m_waveTransmitMat = new Material(Shader.Find("Unlit/WaveTransmitShader"));
        m_waterWaveMarkTexture = new RenderTexture(FFTSize, FFTSize, 0, RenderTextureFormat.Default);
        m_waterWaveMarkTexture.name = "m_waterWaveMarkTexture";
        m_waveTransmitTexture = new RenderTexture(FFTSize, FFTSize, 0, RenderTextureFormat.Default);
        m_waveTransmitTexture.name = "m_waveTransmitTexture";
        m_prevWaveMarkTexture = new RenderTexture(FFTSize, FFTSize, 0, RenderTextureFormat.Default);
        m_prevWaveMarkTexture.name = "m_prevWaveMarkTexture";

        InitWaveTransmitParams();

        TempRT = new RenderTexture(FFTSize, FFTSize, 0, RenderTextureFormat.Default);
        TempRT.Create();
        particlePosTex = new Texture2D(FFTSize, FFTSize, TextureFormat.RFloat, false, true);
        pixData = particlePosTex.GetRawTextureData<float>();
        _waveParticles = WaveParticleSystem.Instance._waveParticles;
        waveFilterMaterial = new Material(waveFilterShader);
        m_TmpHeightFieldRT = CreateRT(FFTSize, "tmp");
        waveParticlePointRT = CreateRT(FFTSize, "particle");

        a.value = Amplitude;
        t.value = TimeScale;
        w.value = WindScale;
    }
    
    private void OceanComputerShaderInit()
    {
        // ocean size, Lx = Lz
        FFTSize = (int)Mathf.Pow(2, FFTPow);

        // create rendering texture
        GaussianRandomRT = CreateRT(FFTSize, "GaussianRandom"); 
        HeightSpectrumRT = CreateRT(FFTSize, "HeightSpectrum");
        DisplaceXSpectrumRT = CreateRT(FFTSize, "DisplaceXSpectrum");
        DisplaceZSpectrumRT = CreateRT(FFTSize, "DisplaceZSpectrum");
        DisplaceRT = CreateRT(FFTSize, "Displace");
        OutputRT = CreateRT(FFTSize, "Output");
        NormalRT = CreateRT(FFTSize, "Normal");
        BubblesRT = CreateRT(FFTSize, "Bubbles");

        // get all kernelID
        kernelComputeGaussianRandom = OceanCS.FindKernel("ComputeGaussianRandom");
        kernelCreateHeightSpectrum = OceanCS.FindKernel("CreateHeightSpectrum");
        kernelCreateDisplaceSpectrum = OceanCS.FindKernel("CreateDisplaceSpectrum");
        kernelFFTHorizontal = OceanCS.FindKernel("FFTHorizontal");
        kernelFFTHorizontalEnd = OceanCS.FindKernel("FFTHorizontalEnd");
        kernelFFTVertical = OceanCS.FindKernel("FFTVertical");
        kernelFFTVerticalEnd = OceanCS.FindKernel("FFTVerticalEnd");
        kernelTextureGenerationDisplace = OceanCS.FindKernel("TextureGenerationDisplace");
        kernelTextureGenerationNormalBubbles = OceanCS.FindKernel("TextureGenerationNormalBubbles");
        // set ComputerShader values
        OceanCS.SetInt("N", FFTSize);
        OceanCS.SetFloat("OceanLength", MeshLength);
        // generate Gaussian random number
        OceanCS.SetTexture(kernelComputeGaussianRandom, "GaussianRandomRT", GaussianRandomRT);
        OceanCS.Dispatch(kernelComputeGaussianRandom, FFTSize / 8, FFTSize / 8, 1);
    }
    
    //create textures
    private RenderTexture CreateRT(int size, string name)
    {
        RenderTexture rt = new RenderTexture(size, size, 0, RenderTextureFormat.ARGBFloat);
        rt.enableRandomWrite = true;
        rt.name = name;
        rt.Create();
        return rt;
    }
    
    private void Update()
    {
        time += Time.deltaTime * TimeScale;
        // compute ocean wave
        ComputeOcean();
        // set texture
        SetMaterialTex();

        WaterMark();
        WaveTransmit();
        //
        generatePosTex();

        Amplitude = a.value;
        TimeScale = t.value;
        WindScale = w.value;

        na.text = $"{Amplitude}";
        nt.text = $"{TimeScale}";
        nw.text = $"{WindScale}";
    }
    
    void generatePosTex()
    {
        if(_waveParticles == null)
            _waveParticles = WaveParticleSystem.Instance._waveParticles;
        if (_waveParticles.Count == 0)
            return;
        for( int i = 0; i < pixData.Length; i++ )
        {
            pixData[ i ] = 0;
        }
        
        float texelW = 1.0f / FFTSize;
        float texelH = 1.0f / FFTSize;
        for (int i = 0; i < _waveParticles.Count; i++)
        {
            Vector3 pos = this.transform.worldToLocalMatrix * _waveParticles[i].data.pos;
            Vector2 posInPlane = new Vector2(pos.x / this.GetComponent<MeshFilter>().mesh.bounds.size.x + 0.5f,
                pos.z / this.GetComponent<MeshFilter>().mesh.bounds.size.z + 0.5f);
            if (posInPlane.x <= 0.01 || posInPlane.y <= 0.01 || posInPlane.x >= 0.99 || posInPlane.y >= 0.99)
                continue;
            // Pixel coordinates with fractional parts
            float xF = posInPlane.x / texelW;
            float yF = posInPlane.y / texelH;
            // Texture pixel indices
            int x = (int)xF;
            int y = (int)yF;
            // Interpolation coefficients between texture indices
            float dX = xF - x;
            float dY = yF - y;
            // Indices 
            int x0y0 = x         + y         * FFTSize;
            int x1y0 = ( x + 1 ) + y         * FFTSize;
            int x0y1 = x         + ( y + 1 ) * FFTSize;
            int x1y1 = ( x + 1 ) + ( y + 1 ) * FFTSize;
            pixData[ x0y0 ] += _waveParticles[i].data.amplitude * ( 1.0f - dX ) * ( 1.0f - dY );
            pixData[ x1y0 ] += _waveParticles[i].data.amplitude * dX            * ( 1.0f - dY );
            pixData[ x0y1 ] += _waveParticles[i].data.amplitude * ( 1.0f - dX ) * dY;
            pixData[ x1y1 ] += _waveParticles[i].data.amplitude * dX            * dY;
        }
        particlePosTex.Apply();
        waveFilterMaterial.SetFloat( "_WaveParticleRadius", _waveParticles[0].data.radius);
        Graphics.Blit( particlePosTex, m_TmpHeightFieldRT, waveFilterMaterial, pass: 0 );
        Graphics.Blit( m_TmpHeightFieldRT, waveParticlePointRT, waveFilterMaterial, pass: 1 ); 
        
        texMaterial.SetTexture("_MainTex", particlePosTex);

    }
    
    //define ocean wave parameters
    private void ComputeOcean()
    {
        // amplitude
        OceanCS.SetFloat("Amplitude", Amplitude);
        // wind speed
        Vector2 wind = new Vector2(Wind.x, Wind.y);
        wind.Normalize();
        wind *= WindScale;
        OceanCS.SetVector("Wind", new Vector2(wind.x, wind.y));
        // set parameters
        OceanCS.SetFloat("Lambda", Lambda);
        OceanCS.SetFloat("HeightScale", HeightScale);
        OceanCS.SetFloat("BubblesScale", BubblesScale);
        OceanCS.SetFloat("BubblesThreshold",BubblesThreshold);
        // generate height texture
        OceanCS.SetFloat("Time", time);
        OceanCS.SetTexture(kernelCreateHeightSpectrum, "GaussianRandomRT", GaussianRandomRT);
        OceanCS.SetTexture(kernelCreateHeightSpectrum, "HeightSpectrumRT", HeightSpectrumRT);
        OceanCS.Dispatch(kernelCreateHeightSpectrum, FFTSize / 8, FFTSize / 8, 1);
        // generate offset texture
        OceanCS.SetTexture(kernelCreateDisplaceSpectrum, "HeightSpectrumRT", HeightSpectrumRT);
        OceanCS.SetTexture(kernelCreateDisplaceSpectrum, "DisplaceXSpectrumRT", DisplaceXSpectrumRT);
        OceanCS.SetTexture(kernelCreateDisplaceSpectrum, "DisplaceZSpectrumRT", DisplaceZSpectrumRT);
        OceanCS.Dispatch(kernelCreateDisplaceSpectrum, FFTSize / 8, FFTSize / 8, 1);
        // horizontal FFT
        for (int m = 1; m <= FFTPow; m++)
        {
            int ns = (int)Mathf.Pow(2, m - 1);
            OceanCS.SetInt("Ns", ns);
            //process with FFT
            if (m != FFTPow)
            {
                ComputeOceanFFT(kernelFFTHorizontal, ref HeightSpectrumRT);
                ComputeOceanFFT(kernelFFTHorizontal, ref DisplaceXSpectrumRT);
                ComputeOceanFFT(kernelFFTHorizontal, ref DisplaceZSpectrumRT);
            }
            else
            {
                ComputeOceanFFT(kernelFFTHorizontalEnd, ref HeightSpectrumRT);
                ComputeOceanFFT(kernelFFTHorizontalEnd, ref DisplaceXSpectrumRT);
                ComputeOceanFFT(kernelFFTHorizontalEnd, ref DisplaceZSpectrumRT);
            }

        }
        // vertical FFT
        for (int m = 1; m <= FFTPow; m++)
        {
            int ns = (int)Mathf.Pow(2, m - 1);
            OceanCS.SetInt("Ns", ns);
            // process with FFT
            if (m != FFTPow)
            {
                ComputeOceanFFT(kernelFFTVertical, ref HeightSpectrumRT);
                ComputeOceanFFT(kernelFFTVertical, ref DisplaceXSpectrumRT);
                ComputeOceanFFT(kernelFFTVertical, ref DisplaceZSpectrumRT);
            }
            else
            {
                ComputeOceanFFT(kernelFFTVerticalEnd, ref HeightSpectrumRT);
                ComputeOceanFFT(kernelFFTVerticalEnd, ref DisplaceXSpectrumRT);
                ComputeOceanFFT(kernelFFTVerticalEnd, ref DisplaceZSpectrumRT);
            }
    
        }
        // offset texture generation 
        OceanCS.SetTexture(kernelTextureGenerationDisplace, "HeightSpectrumRT", HeightSpectrumRT);
        OceanCS.SetTexture(kernelTextureGenerationDisplace, "DisplaceXSpectrumRT", DisplaceXSpectrumRT);
        OceanCS.SetTexture(kernelTextureGenerationDisplace, "DisplaceZSpectrumRT", DisplaceZSpectrumRT);
        OceanCS.SetTexture(kernelTextureGenerationDisplace, "DisplaceRT", DisplaceRT);
        OceanCS.Dispatch(kernelTextureGenerationDisplace, FFTSize / 8, FFTSize / 8, 1);
        // normal texture and foam texture
        OceanCS.SetTexture(kernelTextureGenerationNormalBubbles, "DisplaceRT", DisplaceRT);
        OceanCS.SetTexture(kernelTextureGenerationNormalBubbles, "NormalRT", NormalRT);
        OceanCS.SetTexture(kernelTextureGenerationNormalBubbles, "BubblesRT", BubblesRT);
        OceanCS.Dispatch(kernelTextureGenerationNormalBubbles, FFTSize / 8, FFTSize / 8, 1);
    }

    // compute height of the ocean wave using FFT
    /// <param name="kernel"></param>
    /// <param name="input"></param>
    private void ComputeOceanFFT(int kernel, ref RenderTexture input)
    {
        OceanCS.SetTexture(kernel, "InputRT", input);
        OceanCS.SetTexture(kernel, "OutputRT", OutputRT);
        OceanCS.Dispatch(kernel, FFTSize / 8, FFTSize / 8, 1);
        // exchang input and output texture of the ocean
        RenderTexture rt = input;
        input = OutputRT;
        OutputRT = rt;
    }
    
    //set texture for each material
    private void SetMaterialTex()
    {
        // ocean 
        OceanMaterial.SetTexture("_Displace", DisplaceRT);
        OceanMaterial.SetTexture("_Normal", NormalRT);
        OceanMaterial.SetTexture("_Bubbles", BubblesRT);
        OceanMaterial.SetTexture("_WaveResult", m_waveTransmitTexture);
        OceanMaterial.SetTexture("_Tex", objectRenderTexture);
        // display
        DisplaceXMat.SetTexture("_MainTex", DisplaceXSpectrumRT);
        DisplaceYMat.SetTexture("_MainTex", HeightSpectrumRT);
        DisplaceZMat.SetTexture("_MainTex", DisplaceZSpectrumRT);
        DisplaceMat.SetTexture("_MainTex", DisplaceRT);
        NormalMat.SetTexture("_MainTex", NormalRT);
        BubblesMat.SetTexture("_MainTex", BubblesRT);
        GuassianMat.SetTexture("_MainTex", GaussianRandomRT);
        WaveMat.SetTexture("_MainTex", m_waveTransmitTexture);
    }

    private bool hasHit = false;
    Vector2 hitPos;
    private Vector4 m_waveTransmitParams;
    
    [Header(("interactive Wave"))]
    public float WaveRadius = 0.01f;
    public float WaveSpeed = 1.0f;
    public float WaveViscosity = 1.0f;
    public float WaveAtten = 0.99f; // attenuation
    public float WaveHeight = 0.999f;

    public static OceanVisulizer Instance
    {
        get
        {
            if (sInstance == null)
                sInstance = FindObjectOfType<OceanVisulizer>();
            return sInstance;
        }
    }

    private static OceanVisulizer sInstance;
    
    //collision check
    public void SphereTest(GameObject sphere)
    {
        if (sphere.transform.position.y < 10)
        {
            Vector3 waterPlaneSpacePos = this.transform.worldToLocalMatrix * new Vector4(sphere.transform.position.x, sphere.transform.position.y, sphere.transform.position.z, 1);
            float dx = (waterPlaneSpacePos.x / MeshLength) + 0.5f;
            float dy = (waterPlaneSpacePos.z / MeshLength) + 0.5f;

            m_waveMarkParams.Set(0, 0, WaveRadius * WaveRadius, WaveHeight);
            hasHit = true;
            
            addMaterial.SetTexture("_MainTex", m_waterWaveMarkTexture);
            addMaterial.SetTexture("_Tex", objectRenderTexture);
            Graphics.Blit(null, TempRT, addMaterial);
            RenderTexture rt = TempRT;
            TempRT = m_waterWaveMarkTexture;
            m_waterWaveMarkTexture = rt;
            WaveTransmit();
        }        
        else
        {
            hasHit = false;
        }
    }

    // release space
    /// <exception cref="NotImplementedException"></exception>
    private void OnDestroy()
    {
        GaussianRandomRT.Release();
        HeightSpectrumRT.Release();
        DisplaceXSpectrumRT.Release();
        DisplaceZSpectrumRT.Release();
        DisplaceRT.Release();
        OutputRT.Release();
        NormalRT.Release();
        BubblesRT.Release();
    }
    
    void WaterMark()
    {
        if (hasHit)
        {
            m_waterWaveMarkMat.SetVector("_WaveMarkParams", m_waveMarkParams);
            Graphics.Blit(m_waveTransmitTexture, m_waterWaveMarkTexture, m_waterWaveMarkMat);
        }
    }
    
    void WaveTransmit()
    {
        m_waveTransmitMat.SetVector("_WaveTransmitParams", m_waveTransmitParams);
        m_waveTransmitMat.SetFloat("_WaveAtten", WaveAtten);
        m_waveTransmitMat.SetTexture("_PrevWaveMarkTex", m_prevWaveMarkTexture);

        RenderTexture rt = RenderTexture.GetTemporary(FFTSize, FFTSize);
        Graphics.Blit(m_waterWaveMarkTexture, rt, m_waveTransmitMat);
        Graphics.Blit(m_waterWaveMarkTexture, m_prevWaveMarkTexture);
        Graphics.Blit(rt, m_waterWaveMarkTexture);
        Graphics.Blit(rt, m_waveTransmitTexture);
        RenderTexture.ReleaseTemporary(rt);
    }

    void InitWaveTransmitParams()
    {
        float uvStep = 1.0f / FFTSize;
        float dt = Time.fixedDeltaTime;
        // maximum progressive visosity
        float maxWaveStepVisosity = uvStep / (2 * dt) * (Mathf.Sqrt(WaveViscosity * dt + 2));
        // visosity squared u^2
        float waveVisositySqr = WaveViscosity * WaveViscosity;
        // current speed
        float curWaveSpeed = maxWaveStepVisosity * WaveSpeed;
        // speed squared c^2
        float curWaveSpeedSqr = curWaveSpeed * curWaveSpeed;
        // single wave transform squared d^2
        float uvStepSqr = uvStep * uvStep;

        float i = Mathf.Sqrt(waveVisositySqr + 32 * curWaveSpeedSqr / uvStepSqr);
        float j = 8 * curWaveSpeedSqr / uvStepSqr;

        // wave transfer equation
        // (4 - 8 * c^2 * t^2 / d^2) / (u * t + 2) + (u * t - 2) / (u * t + 2) * z(x,y,z, t - dt) + (2 * c^2 * t^2 / d ^2) / (u * t + 2)
        // * (z(x + dx,y,t) + z(x - dx, y, t) + z(x,y + dy, t) + z(x, y - dy, t);

        //ut
        float ut = WaveViscosity * dt;
        //c^2 * t^2 / d^2
        float ctdSqr = curWaveSpeedSqr * dt * dt / uvStepSqr;
        // ut + 2
        float utp2 = ut + 2;
        // ut - 2
        float utm2 = ut - 2;
        //(4 - 8 * c^2 * t^2 / d^2) / (u * t + 2) 
        float p1 = (4 - 8 * ctdSqr) / utp2;
        //(u * t - 2) / (u * t + 2)
        float p2 = utm2 / utp2;
        //(2 * c^2 * t^2 / d ^2) / (u * t + 2)
        float p3 = (2 * ctdSqr) / utp2;

        m_waveTransmitParams.Set(p1, p2, p3, uvStep);

        //Debug.LogFormat("i {0} j {1} maxSpeed {2}", i, j, maxWaveStepVisosity);
        //Debug.LogFormat("p1 {0} p2 {1} p3 {2}", p1, p2, p3);
    }
    
}